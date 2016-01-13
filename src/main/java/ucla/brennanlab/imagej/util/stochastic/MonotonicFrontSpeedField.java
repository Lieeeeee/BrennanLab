package ucla.brennanlab.imagej.util.stochastic;

import cern.jet.random.Normal;
import ij.IJ;
import org.ujmp.core.Matrix;
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
    private double[][] currentMean;
    private double[][] currentVariance;

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
        this.currentMean = new double[width][height];
        this.currentVariance = new double[width][height];
        for(int i=0; i<width;i++){
            for(int j=0;j<height;j++){
                currentMean[i][j] = priormean;
                currentVariance[i][j] = priorvar;
            }
        }

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
        double totspeed = 0;
        double totspeed2 = 0;
        for(int j=0; j<incomingBoundaryCoordinates.length;j++){
            incomingBoundarySpeeds[j] = (float) Math.abs(
                    this.wavePositions
                            .get(this.wavePositions.size()-1)
                            .get(incomingBoundaryCoordinates[j][0],incomingBoundaryCoordinates[j][1])/
                    (time-times.get(this.wavePositions.size()-1)));
            totspeed += incomingBoundarySpeeds[j];
            totspeed2 += Math.pow(incomingBoundarySpeeds[j],2);
        }
        this.currentMean = new double[width][height];
        this.currentVariance = new double[width][height];
        for(int i=0; i<width;i++){
            for(int j=0;j<height;j++){
                currentMean[i][j] = totspeed/incomingBoundaryCoordinates.length;
                currentVariance[i][j] = totspeed2/incomingBoundaryCoordinates.length-Math.pow(currentMean[i][j],2);
            }
        }
        IJ.log("boundary speed: "+totspeed/incomingBoundaryCoordinates.length);

        double[][] covariates = new double[1][1];
        covariates[0][0] = 1;
        //krigingLattice.addObservations(incomingBoundaryCoordinates,covariates,incomingBoundarySpeeds);
        this.times.add(time);
        this.wavePositions.add(incomingShape);

        //computeCurrentMean();
        //computeCurrentVariance();

    }

 public MonotonicFrontSpeedField clone() {
        MonotonicFrontSpeedField s = new MonotonicFrontSpeedField(width, height);
        s.krigingLattice = this.krigingLattice;
        s.priorVarSpeed = priorVarSpeed;
        s.priorMeanSpeed = priorMeanSpeed;

        return s;
    }


    /**
     * Sample a speed field
     * @return
     */
    public double[][] sample(ImplicitShape2D reference, int band){
        /**
         * Make a list of the coordinates
         */

        /**
         * Assemble V0
         */

        return null;
    }

    public ArrayList<double[][]> sample(int samples,ImplicitShape2D reference, int band){
        ArrayList<double[][]> out = new ArrayList<double[][]>(samples);
        for(int j=0;j<samples;j++){
            out.add(sample(reference,band));
        }
        return out;
    }

    /**
     * Return the current mean estimate for the speed field
     * @return
     */
    public double[][] computeCurrentMean() {
        /**
         * Retrieve \hat\beta from kriging object and fill out missing locations
         * with the mean speed
         */
        double[][] s = new double[width][height];

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
            Matrix beta = krigingLattice.getBeta();

        }

        return s;
    }

    /**
     *
     * @return
     */
    public double[][] getCurrentMean(){
        return this.currentMean;
    }

    public double[][] getCurrentVariance(){
        return this.currentVariance;
    }

    /**
     * Return the current variance estimate for the speed field
     * @return
     */
    public double[][] computeCurrentVariance() {
        /**
         * Retrieve \hat\beta from kriging object and fill out missing locations
         * with the mean speed
         */

        double[][] s = new double[width][height];

        /**
         * If kriging object is null, we do not have any observations yet, so
         * draw according to the prior
         */
        if (krigingLattice == null) {

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    s[x][y] = (float) this.priorVarSpeed;
                }
            }

        } else {
            Matrix beta = krigingLattice.getBetaVariance();

        }
        return s;
    }





}
