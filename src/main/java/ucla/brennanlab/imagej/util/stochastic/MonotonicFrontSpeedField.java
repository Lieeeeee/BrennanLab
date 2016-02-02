package ucla.brennanlab.imagej.util.stochastic;

import cern.jet.random.Normal;
import org.delaunay.algorithm.Triangulation;
import org.ujmp.core.Matrix;
import ucla.brennanlab.DelaunayInterpolator;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;

import java.util.ArrayList;


/**
 * This class stores interfaces with the Kriging class to interpolate arrival times
 * and sample speed fields. The kringing lattice is computed on a sublattice
 * <p>
 * Fundamentally, we probabilistically solve the eikonal equation |\nabla T| = F, where $T$
 * are the arrival times.
 *
 * @author Josh Chang
 */
public class MonotonicFrontSpeedField {
    public int width, height;
    public ArrayList<ImplicitShape2D> wavePositions;
    ArrayList<Float> times;
    Normal standardNormal; // use this to generate iid Gaussian variables
    private Kriging2DLattice krigingLattice; // stores actual speeds that we will use
    private double priorMeanSpeed;
    private double priorVarSpeed;
    public double overallMeanSpeed; // Temporary
    public double overallSDSpeed; // Temporary
    private double[][] currentMean;
    private double[][] currentVariance;
    private int xreduction;
    private int yreduction;

    /**
     * @param width  Width of interpolation signedDistance
     * @param height Width of height signedDistance
     **/
    public MonotonicFrontSpeedField(int width, int height) {
        this(width, height, Float.MIN_VALUE, Float.MAX_VALUE, 1, 1);
    }

    public MonotonicFrontSpeedField(int width, int height, double priormean, double priorvar ,int xreduction, int yreduction) {
        this.width = width;
        this.height = height;
        this.times = new ArrayList<Float>();
        this.wavePositions = new ArrayList<ImplicitShape2D>();
        this.priorMeanSpeed = priormean;
        this.priorVarSpeed = priorvar;
        this.standardNormal = new Normal(0, 1,
                new cern.jet.random.engine.MersenneTwister(new java.util.Date())
        );
        krigingLattice = new Kriging2DLattice(priormean,1.0/priorvar);
        this.currentMean = new double[width][height];
        this.currentVariance = new double[width][height];
        for(int i=0; i<width;i++){
            for(int j=0;j<height;j++){
                currentMean[i][j] = priormean;
                currentVariance[i][j] = priorvar;
            }
        }
        this.xreduction = xreduction;
        this.yreduction = yreduction;

    }


    /**
     * Add an arrival by mask
     * @param mask
     * @param time
     */
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
            incomingBoundarySpeeds[j] = (float) Math.sqrt( Math.abs(
                    this.wavePositions
                            .get(this.wavePositions.size()-1)
                            .get(incomingBoundaryCoordinates[j][0],incomingBoundaryCoordinates[j][1])/
                    (time-times.get(this.wavePositions.size()-1))));
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

        double[][] covariates = new double[incomingBoundaryCoordinates.length][1];
        for(int j=0;j<incomingBoundaryCoordinates.length;j++){
            covariates[j][0] = 1;
        }
        //addObservations(incomingBoundaryCoordinates,covariates,incomingBoundarySpeeds);
        this.times.add(time);
        this.wavePositions.add(incomingShape);

        this.overallMeanSpeed = totspeed/incomingBoundaryCoordinates.length;
        this.overallSDSpeed = Math.sqrt(totspeed2/incomingBoundaryCoordinates.length-Math.pow(this.overallMeanSpeed,2));


    }

    /***
     * Ignore this. stuff and just infer using only the inforamtion given
     * @param locations
     * @param covariates
     * @param speeds
     * @return
     */
    public float[][] locallyInferAndSample(int[][] locations, double[][] covariates, double[] speeds){

        Kriging2DLattice kriger = new Kriging2DLattice(this.priorMeanSpeed,1.0/this.priorVarSpeed);
        kriger.addObservations(locations,covariates,speeds);

        return null;

    }

    public MonotonicFrontSpeedField clone() {
        MonotonicFrontSpeedField s = new MonotonicFrontSpeedField(width, height);
        s.krigingLattice = this.krigingLattice;
        s.priorVarSpeed = priorVarSpeed;
        s.priorMeanSpeed = priorMeanSpeed;

        return s;
    }

    public void addObservations(int[][] locations, double[][] covariates, double[] speeds){
        int[][] incomingBoundaryCoordinates = new int[locations.length][locations[0].length];



        for(int i=0; i<locations.length;i++){
            incomingBoundaryCoordinates[i][0] = getXcoord(locations[i][0]);
            incomingBoundaryCoordinates[i][1] = getYcoord(locations[i][1]);
        }
        krigingLattice.addObservations(incomingBoundaryCoordinates,covariates,speeds);

    }

    public ArrayList<double[][]> sampleUnscented(double alpha, double beta, double lambda, int L){
        ArrayList<double[][]> samples = new ArrayList<double[][]>(2*L+1);
        double[][] mean = computeCurrentMeanField();
        double[][] variance = computeCurrentVarianceField();


        for(int j=1; j<= L; j++){

        }
        return samples;
    }

    private double[][] matrixPlus(double[][] a, double[][] b){
        double[][] out = new double[a.length][a[0].length];
        for(int i=0;i<a.length;i++){
            for(int j=0;j<a[0].length;j++){
                out[i][j] = a[i][j] + b[i][j];
            }
        }
        return out;
    }

    private double[][] matrixSqrt(double[][] a){
        double[][] out = new double[a.length][a[0].length];
        for(int i=0;i<a.length;i++){
            for(int j=0;j<a[0].length;j++){
                out[i][j] = Math.sqrt(a[i][j]);
            }
        }
        return out;
    }

    /**
     * Return the current mean estimate for the speed field
     * @return
     */
    public double[][] computeCurrentMeanField() {
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

    public double[][] getCurrentVarianceField(){
        return this.currentVariance;
    }


    /**
     * Return the current variance estimate for the speed field
     * @return
     */
    public double[][] computeCurrentVarianceField() {
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



    public float[][] interpolateSpeedField(){
        float[][] speedgrid = new float[width][height];


        return null;
    }

    public float[][] interpolateSpeedFromLevelSets(
            ArrayList<ImplicitShape2D> levelsetSignedDistances) {
        DelaunayInterpolator t = new DelaunayInterpolator();
        float[][] speedgrid;

        for(int j=1; j<levelsetSignedDistances.size();j++){
            int[][] boundaryCoordinates = levelsetSignedDistances.get(j).getBoundaryCoordinates();
            float[] boundarySpeeds = new float[boundaryCoordinates.length];
            for(int k=0; k<boundaryCoordinates.length;k++){
                int x = boundaryCoordinates[k][0];
                int y = boundaryCoordinates[k][1];
                boundarySpeeds[k] = (float)
                        Math.sqrt(Math.abs((levelsetSignedDistances.get(j-1).get(x,y))));
                t.addVertex((double) x,(double) y,(double) boundarySpeeds[k]);
            }
        }
        try {
            t.triangulate();
        } catch (Triangulation.InvalidVertexException e) {
            e.printStackTrace();
        }
        speedgrid = t.getInterpolation(width,height);
        return speedgrid;
    }

    public float[][] listBoundarySpeedsFromLevelSets(ArrayList<ImplicitShape2D> levelsetSignedDistances){
        ArrayList<Integer> xcoords = new ArrayList<Integer>();
        ArrayList<Integer> ycoords = new ArrayList<Integer>();
        ArrayList<Double> speed = new ArrayList<Double>();
        ArrayList<Double> times = new ArrayList<Double>();
        for(int j=1; j<levelsetSignedDistances.size();j++){
            int[][] boundaryCoordinates = levelsetSignedDistances.get(j).getBoundaryCoordinates();
            float[] boundarySpeeds = new float[boundaryCoordinates.length];
            for(int k=0; k<boundaryCoordinates.length;k++){
                int x = boundaryCoordinates[k][0];
                xcoords.add(x);
                int y = boundaryCoordinates[k][1];
                ycoords.add(y);
                boundarySpeeds[k] = (float)
                        Math.sqrt(Math.abs((levelsetSignedDistances.get(j-1).get(x,y))));
                speed.add((double)boundarySpeeds[k]);
                times.add((double)this.times.get(j));
            }
        }
        float[][] out = new float[xcoords.size()][4];
        for(int i=0;i<xcoords.size();i++){
            out[i][0] = xcoords.get(i);
            out[i][1] = ycoords.get(i);
            out[i][2] = speed.get(i).floatValue();
            out[i][3] = times.get(i).floatValue();
        }
        return out;

    }



    public int getXcoord(int x){
        return (int) Math.floor(x/this.xreduction);
    }

    public int getYcoord(int y){
        return (int)Math.floor(y/this.yreduction);
    }
}
