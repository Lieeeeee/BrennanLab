package ucla.brennanlab.imagej.util.stochastic;

import cern.jet.random.Normal;
import ij.IJ;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;

import java.util.ArrayList;

/**
 * This class stores interfaces with the Kriging class to interpolate arrival times
 * and sample speed fields.
 * <p>
 * Fundamentally, we probabilistically solve the eikonal equation |\nabla T| = F, where $T$
 * are the arrival times.
 *
 * @author Josh Chang
 */
public class MonotonicFrontSpeedField {
    public int width, height;
    ArrayList<ImplicitShape2D> wavePositions;
    ArrayList<Float> times;
    Normal standardNormal; // use this to generate iid Gaussian variables
    private Kriging2DLattice krigingLattice; // stores actual speeds that we will use
    private double priorMeanSpeed;
    private double priorVarSpeed;

    /**
     * @param width  Width of interpolation signedDistance
     * @param height Width of height signedDistance
     **/
    public MonotonicFrontSpeedField(int width, int height) {
        this(width, height, Float.MIN_VALUE, Float.MAX_VALUE);
    }

    public MonotonicFrontSpeedField(int width, int height, double priormean, double priorvar) {
        this.width = width;
        this.height = height;
        this.times = new ArrayList<Float>();
        this.wavePositions = new ArrayList<ImplicitShape2D>();
        this.priorMeanSpeed = priormean;
        this.priorVarSpeed = priorvar;
        this.standardNormal = new Normal(
                0,
                1,
                new cern.jet.random.engine.MersenneTwister(new java.util.Date()));
        krigingLattice = new Kriging2DLattice(priormean,1.0/priorvar);
    }


    public void addArrival(boolean[][] mask, float time){
        this.addArrival(new ImplicitShape2D(mask), time);


    }

    /**
     * Add another ROI and update the kriging object. We don't perform any
     * sanity checking here.
     *
     * @param incomingShape The arrival times with null for not-defined
     */
    public void addArrival(ImplicitShape2D incomingShape, float time) {

        /**
         * Check to make sure that:
         *      1) New front does not cross the other observed fronts inadmissably
         *      2) ... anything else?
         */

        if(this.times.size()==0){
            this.times.add(time);
            this.wavePositions.add(incomingShape);
            return;
        }

        // Figure out position to insert into this.times

        // For now assume that we are adding times in order!



        /**
         * For each point on the front of incomingShape, calculate the signed
         * distance from the previous shape. Use this computation to approximate the
         * speed
         */

        int[][] incomingBoundaryCoordinates = incomingShape.getBoundaryCoordinates();
        double[] incomingBoundarySpeeds = new double[incomingBoundaryCoordinates.length];
        for(int j=0; j<incomingBoundaryCoordinates.length;j++){
            incomingBoundarySpeeds[j] = (float) Math.abs(
                    this.wavePositions
                            .get(this.wavePositions.size())
                            .get(incomingBoundaryCoordinates[j][0],incomingBoundaryCoordinates[j][1])/
                    (time-times.get(this.wavePositions.size())));
        }

        double[][] covariates = new double[1][1];
        covariates[0][0] = 1;
        krigingLattice.addObservations(incomingBoundaryCoordinates,covariates,incomingBoundarySpeeds);
        this.times.add(time);
        this.wavePositions.add(incomingShape);

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

        IJ.log("Mean linearly interpolated speed was "+ tot/n);
        return speeds;
    }

    public MonotonicFrontSpeedField clone() {
        MonotonicFrontSpeedField s = new MonotonicFrontSpeedField(width, height);
        s.setKrigingObject(this.krigingLattice);
        s.priorVarSpeed = priorVarSpeed;
        s.priorMeanSpeed = priorMeanSpeed;

        return s;
    }


    /**
     *
     */

    public Kriging2DLattice getKrigingObject() {
        return krigingLattice;
    }

    public void setKrigingObject(Kriging2DLattice krig) {
        this.krigingLattice = krig;
    }

    public float[][] getPosteriorMeanSpeedField() {
        // poll the kriging object to get the mean

        float[][] s = new float[width][height];

        /**
         * If kriging object is null, we do not have any observations yet, so
         * draw according to the prior
         */
        if (krigingLattice == null) {

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    s[x][y] = (float) this.priorMeanSpeed;
                }
            }

        } else {

        }

        return s;
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
        double[] center = levelsetSignedDistances.get(0).inertialCenter;
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


    /**
     * Sample a speed field
     * @return
     */
    public double[][] sample(){


        return null;
    }

    /**
     * Return the current mean estimate for the speed field
     * @return
     */
    public double[][] currentMean() {
        /**
         * Retrieve \hat\beta from kriging object and fill out missing locations
         * with the mean speed
         */

        // Store as this.runningMean
        return null;
    }

    /**
     * Return the current variance estimate for the speed field
     * @return
     */
    public double[][] currentVariation() {
        /**
         * Retrieve \hat\beta from kriging object and fill out missing locations
         * with the mean speed
         */

        // Store as this.runningVariation
        return null;
    }


    public void updateGMRF() {
        // TODO Auto-generated method stub

    }

    /**
     * Return the interpolated arrival times implied by the probability model
     * @return
     */
    public float[][] getArrivalTimes(){
        return null;
    }



}
