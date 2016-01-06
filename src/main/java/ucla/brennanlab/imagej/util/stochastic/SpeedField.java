package ucla.brennanlab.imagej.util.stochastic;

import cern.jet.random.Normal;
import ij.IJ;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;

import java.util.ArrayList;

/**
 * This class stores a single realization of speed field, and includes methods
 * for interpolation. The speed field calls upon Kriging2D field to interpolate
 * the speeds in a signedDistance. "Regression to the mean" is expected, so
 * speeds far away from observed points are simply assigned mean-speeds for ease
 * of computation (To save time)
 * <p>
 * Also stores the kriging object and does translation between it and image
 * spatial field. Of note is the fact that the representation of speed here is
 * unitless, it is expected that there is 1 speed unit between each ROI passed
 * in
 *
 * @author Joshua Chang {@link mailto:joshchang@ucla.edu}
 */
public class SpeedField {
    public int width, height;
    public double latestSpeed;
    ArrayList<ImplicitShape2D> wavePositions;
    Normal standardNormal;
    private float[][] speedField; // "response" associated with wavePositions
    private float[][] arrivalTimes; // arrival times of wavePositions
    private Kriging2DLattice krig; // stores actual speeds that we will use
    private float priorMeanSpeed;
    private float priorVarSpeed;

    /**
     * @param width  Width of interpolation signedDistance
     * @param height Width of height signedDistance
     **/
    public SpeedField(int width, int height) {
        this(width, height, Float.MIN_VALUE, Float.MAX_VALUE);
    }

    public SpeedField(int width, int height, float priormean, float priorvar) {
        this.width = width;
        this.height = height;
        this.speedField = new float[width][height];
        this.wavePositions = new ArrayList<ImplicitShape2D>();
        this.priorMeanSpeed = priormean;
        this.setPriorVarSpeed(priorvar);
        this.standardNormal = new Normal(
                0,
                1,
                new cern.jet.random.engine.MersenneTwister(new java.util.Date()));
    }

    /**
     * Add another ROI and update the kriging object. We don't perform any
     * sanity checking here.
     *
     * @param thisLevelset The arrival times with null for not-defined
     */
    public void addArrival(ImplicitShape2D thisLevelset) {

        if (wavePositions.size() == 0) {
            wavePositions.add(thisLevelset);
            /**
             * This is the first arrival, so just find the center and add it
             */
            this.arrivalTimes = interpolateArrivalsFromLevelsets(wavePositions);
            this.speedField = arrivalsToSpeed(this.arrivalTimes);

        } else {
            ImplicitShape2D prevLevelset = wavePositions.get(wavePositions
                    .size() - 1);

            double totTime = 0;
            int n = 0;
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    if (thisLevelset.get(x, y) <= 0
                            && prevLevelset.get(x, y) > 0) {
                        // interpolate the arrival times
                        this.arrivalTimes[x][y] = (float) ((-this.wavePositions.size()
                                * thisLevelset.get(x, y) + (this.wavePositions
                                .size() + 1)
                                * prevLevelset.get(x, y))
                                / (prevLevelset.get(x, y) - thisLevelset.get(x,
                                y)));
                        n++;
                        totTime += this.arrivalTimes[x][y];
                    }
                }
            }

            IJ.log("Average interpolated arrival time: " + totTime / n);

            double totspeed = 0;
            n = 0;
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    if (thisLevelset.get(x, y) <= 0
                            && prevLevelset.get(x, y) > 0) {
                        // interpolate the upwind finite differences speed!
                        int upX = Math.min(x + 1, width - 1);
                        int upY = Math.min(y + 1, height - 1);
                        int downX = Math.max(x - 1, 0);
                        int downY = Math.max(y - 1, 0);
                        double dTx = Math.max(Math.max(this.arrivalTimes[x][y] -

                                this.arrivalTimes[upX][y], this.arrivalTimes[x][y]
                                - this.arrivalTimes[downX][y]), 0);
                        double dTy = Math.max(Math.max(this.arrivalTimes[x][y]
                                        - this.arrivalTimes[x][upY],
                                this.arrivalTimes[x][y]
                                        - this.arrivalTimes[x][downY]), 0);

                        // Maybe use second order scheme?

                        this.speedField[x][y] = (float) Math.pow(dTx * dTx
                                + dTy * dTy, -0.5);
                        if (x > 5 && y > 5 && x < width - 5
                                && y < width - 5) {
                            totspeed += speedField[x][y];
                            n++;
                        }
                    }
                }
            }
            if (!Double.isInfinite(totspeed) && !(n == 0)) {
                this.latestSpeed = totspeed / n;
            }

            wavePositions.add(thisLevelset);
            IJ.log("Average wavespeed " + totspeed / n + " pixels/frame");
        }

    }

    /**
     * Norm the gradient of the arrival times, then invert
     *
     * @param arrivals Matrix of arrival times
     * @return
     */
    public float[][] arrivalsToSpeed(float[][] arrivals) {
        int n = 0;
        float tot = 0;
        float[][] speeds = new float[width][height];
        float dTx, dTy;

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (Float.isInfinite(arrivals[x][y])
                        || Float.isNaN(arrivals[x][y]))
                    speeds[x][y] = Float.NaN;
                else {

                    int upX = Math.min(x + 1, width - 1);
                    int upY = Math.min(y + 1, height - 1);
                    int downX = Math.max(x - 1, 0);
                    int downY = Math.max(y - 1, 0);
            /*
		     * dTx = Math.max(Math.max(arrivals[x][y] -
		     * arrivals[upX][y], arrivals[x][y] - arrivals[downX][y]),
		     * 0); dTy = Math.max(Math.max(arrivals[x][y] -
		     * arrivals[x][upY], arrivals[x][y] - arrivals[x][downY]),
		     * 0);
		     */

                    dTx = (float) (0.5 * (arrivals[upX][y] - arrivals[downX][y]));
                    dTy = (float) (0.5 * (arrivals[x][upY] - arrivals[x][downY]));
                    /*** */

                    /**
                     * Second order finite difference
                     *
                     *
                     * dTx = (float)
                     * (-0.5*(arrivals[upX][y]-4*arrivals[x][y]+3*arrivals
                     * [downX][y])); dTy = (float)
                     * (-0.5*(arrivals[x][upY]-4*arrivals
                     * [x][y]+3*arrivals[x][downY]));/*
                     */

                    speeds[x][y] = (float) Math
                            .pow(dTx * dTx + dTy * dTy, -0.5);
                    if (Float.isInfinite(speeds[x][y])) {
                        speeds[x][y] = Float.NaN;
                    } else if (arrivals[x][y] < 4 && arrivals[x][y] > 2) {
                        tot += speeds[x][y];
                        n++;

                    }

                }
            }
        }

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (arrivals[x][y] < 1)
                    speeds[x][y] = tot / n;
            }
        }

        // IJ.showMessage("Mean linearly interpolated speed was "+ tot/n);
        return speeds;
    }

    public SpeedField clone() {
        SpeedField s = new SpeedField(width, height);
        s.setKrigingObject(this.krig);
        s.setspeedField(this.speedField.clone());
        s.setarrivalTimes(this.arrivalTimes.clone());
        s.setPriorVarSpeed(priorVarSpeed);
        s.setPriorMeanSpeed(priorMeanSpeed);

        return s;
    }

    /**
     * Calculates the coordinates of the center of an ROI TODO fix this when
     * origin is probably not in the boundary of the image, in that case this
     * should give a coordinate on the boundary at the very least
     *
     * @param r
     * @return
     */
    public float[] getCenter(ImplicitShape2D signedDistance) {
        if (signedDistance.get(0, 0) > width * height + 1)
            throw new IllegalArgumentException("ROI is null for some reason");
        float center[] = new float[2];
        int tx = 0, ty = 0, n = 0, x = 0, y = 0;
        while (x < width) {
            for (y = 0; y < height; y++) {
                if (signedDistance.get(x, y) < 0) {
                    tx += x;
                    ty += y;
                    n++;
                }

            }
            x++;
        }
        center[0] = tx / n;
        center[1] = ty / n;
        IJ.log("CSD origin at x: " + center[0] + " y: " + center[1]);
        return (center);
    }

    /**
     *
     */

    public Kriging2DLattice getKrigingObject() {
        return krig;
    }

    public void setKrigingObject(Kriging2DLattice krig) {
        this.krig = krig;
    }

    public float[][] getPosteriorMeanSpeedField() {
        // poll the kriging object to get the mean
        // TODO Auto-generated method stub

        float[][] s = new float[width][height];

        /**
         * If kriging object is null, we do not have any observations yet, so
         * draw according to the prior
         */
        if (krig == null) {

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    s[x][y] = this.priorMeanSpeed;
                }
            }

        } else {

        }

        return s;
    }

    public float getPriorMeanSpeed() {
        return priorMeanSpeed;
    }

    public void setPriorMeanSpeed(float priorMeanSpeed) {
        this.priorMeanSpeed = priorMeanSpeed;
    }

    public float getPriorVarSpeed() {
        return priorVarSpeed;
    }

    public void setPriorVarSpeed(float priorVarSpeed) {
        this.priorVarSpeed = priorVarSpeed;
    }

    public float[][] getspeedField() {

        return this.speedField;
    }

    /**
     * Calculates the arrival time T(x,y) from the sequence of ROIs representing
     * position. This function calculates the arrival times between two ROIs as
     * <p>
     * T = (d_2t_1 + d_1 t_2)/(d_1+d_2), where d_1 is distance from contour 1
     * and d_2 is distance from contour 2. This interpolation is probably
     * justifiable??
     *
     * @param levelsetSignedDistances An array of boundaries from which to calculate the arrival
     *                                times
     * @return
     */
    public float[][] interpolateArrivalsFromLevelsets(
            ArrayList<ImplicitShape2D> levelsetSignedDistances) {
        if (levelsetSignedDistances == null
                || levelsetSignedDistances.size() == 0)
            throw new IllegalArgumentException(
                    "Must have at least a single ROI");
        float[][] arrivalGrid = new float[width][height];
        // Compute center of first roi, set that value to 1
        // TODO back-interpolate to find true start time based on speed
        ImplicitShape2D rls, rlsNext;
        float[] center = getCenter(levelsetSignedDistances.get(0));
        rls = levelsetSignedDistances.get(0);

        // find the center
        int x = 0, y = 0;
        for (x = 0; x < width; x++) {
            for (y = 0; y < height; y++) {

                if (rls.get(x, y) <= 0) {
                    float d = (float) Math.abs(rls.get(x, y));
                    float d0 = (float) Math.sqrt((x - center[0])
                            * (x - center[0]) + (y - center[1])
                            * (y - center[1]));
                    arrivalGrid[x][y] = d0 / (d + d0);
                } else
                    arrivalGrid[x][y] = Float.POSITIVE_INFINITY;
            }
        }

        for (int j = 1; j < levelsetSignedDistances.size() - 1; j++) {
            rlsNext = levelsetSignedDistances.get(j);
            for (x = 0; x < width; x++) {
                for (y = 0; y < height; y++) {
                    if (rlsNext.get(x, y) <= 0 & rls.get(x, y) > 0) {
                        float d1 = (float) Math.abs(rls.get(x, y));
                        float d2 = (float) Math.abs(rlsNext.get(x, y));
                        arrivalGrid[x][y] = (d2 * (j + 1) + d1 * (j + 2))
                                / (d1 + d2);
                    }
                }
            }
            rls = rlsNext;
        }
        return arrivalGrid;
    }

    /**
     *
     */
    public void recalculateKrigingObject() {

    }

    public void setarrivalTimes(float[][] a) {
        this.arrivalTimes = a;
    }

    /**
     * Set the speed field
     *
     * @param speedField
     */
    public void setspeedField(float[][] speedField) {
        this.speedField = speedField;

    }

    /**
     *
     */

    public float[][] getArrivalTimes() {
        // TODO Auto-generated method stub
        return this.arrivalTimes;
    }

    public void updateGMRF() {
        // TODO Auto-generated method stub

    }

}
