package org.brennanlab.imagej.util.stochastic;

import cern.jet.random.Normal;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.apache.commons.math3.analysis.interpolation.PiecewiseBicubicSplineInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.PiecewiseBicubicSplineInterpolator;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.brennanlab.imagej.util.FourierConvolution;
import org.brennanlab.imagej.util.levelsets.ImplicitShape2D;

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
    public double betahat; // Temporary
    public double Abeta; // Temporary
    private double[][] currentMean;
    private double[][] coarseMean;
    private int xreduction;
    private int yreduction;

    /**
     * @param width  Width of interpolation signedDistance
     * @param height Width of height signedDistance
     **/
    public MonotonicFrontSpeedField(int width, int height) {
        this(width, height, Float.MIN_VALUE, Float.MAX_VALUE, 1, 1);
    }

    /**
     *
     * @param width      width of field
     * @param height     height of field
     * @param priormean  mean speed guess
     * @param priorvar   variance of guess
     * @param xreduction for reducing the grid in the x direction
     * @param yreduction reduction in the y direction
     */
    public MonotonicFrontSpeedField(int width, int height, double priormean, double priorvar ,int xreduction, int yreduction) {
        this.width = width;
        this.height = height;
        this.times = new ArrayList<Float>();
        this.wavePositions = new ArrayList<ImplicitShape2D>();
        this.priorMeanSpeed = priormean;
        this.priorVarSpeed = priorvar;
        this.betahat = priormean;
        this.standardNormal = new Normal(0, 1,
                new cern.jet.random.engine.MersenneTwister(new java.util.Date())
        );
        krigingLattice = new Kriging2DLattice(priormean,1.0/priorvar,(int)Math.ceil(width/xreduction),(int)Math.ceil(height/yreduction));
        this.coarseMean = new double[(int)Math.ceil(width/xreduction)][(int)Math.ceil(height/yreduction)];
        this.currentMean = new double[width][height];
        for(int y=0;y<height;y++){
            for(int x=0; x<width;x++){
                currentMean[x][y] = priorMeanSpeed;
            }
        }
        this.xreduction = xreduction;
        this.yreduction = yreduction;

    }


    /**
     * Add an arrival by mask
     * @param mask  boolean mask describing the shap
     * @param time  arrival time
     */
    public void addArrival(boolean[][] mask, float time){
        this.addArrival(new ImplicitShape2D(mask), time);
    }

    /**
     * Add another ROI and update the kriging object. We don't perform any
     * sanity checking here.
     *
     * @param incomingShape The shape that has arrived
     * @param time The time these arrivals have arrive
     */
    public void addArrival(ImplicitShape2D incomingShape, float time) {

        /**
         * @TODO Check to make sure that:
         *      1) New front does not cross the other observed fronts inadmissably
         *      2) ... anything else?
         * @TODO Should really compute the speed field self-consistenly from regularization by
         * solving the appropriate Euler-Lagrange equations.
         */

        if(this.times.size()==0){
            this.times.add(time);
            this.wavePositions.add(incomingShape);
            return;
        }

        // Figure out position to insert into this.times
        int position = this.times.size();
        while(time <= this.times.get(position-1) && position>0){
            position--;
        }

        /**
         * For each point on the front of incomingShape, calculate the signed
         * distance from the previous shape. Use this computation to approximate the
         * speed
         */

        int[][] incomingBoundaryCoordinates = incomingShape.getBoundaryCoordinates();

        ArrayList<int[]> coords = new ArrayList<int[]>();
        ArrayList<Double> spd = new ArrayList<Double>();

        for(int j=0; j<incomingBoundaryCoordinates.length;j++){
            int x = incomingBoundaryCoordinates[j][0];
            int y = incomingBoundaryCoordinates[j][1];
            if(x<=1 || y<=1 || x>=width-2 || y>=height-2) continue;
            coords.add(new int[]{x,y});
            spd.add( Math.abs(
                    wavePositions.get(wavePositions.size()-1)
                            .get(x,y)));

        }

        int[][] coords2 = new int[coords.size()][2];
        double[] speeds = new double[coords.size()];
        for(int i=0;i<coords.size();i++){
            coords2[i][0] = coords.get(i)[0];
            coords2[i][1] = coords.get(i)[1];
            speeds[i] = spd.get(i);
        }

        addObservations(coords2,speeds);
        times.add(position,time);
        wavePositions.add(position,incomingShape);

    }

    public void addObservations(int[][] coordinates, double[] speeds){
        // First condense the observations into coarse grid
        // @TODO Use finite differencing interpolation (bilinear interpolation)
        assert(coordinates.length==speeds.length);

        /**
         * aggregate coordinates, grouping by coarsened coordinate
         */

        ArrayList<int[]> existingpoints = krigingLattice.getPoints();
        ArrayList<Integer> counts = new ArrayList<Integer>();
        ArrayList<int[]> aggregatedpoints = new ArrayList<int[]>();
        ArrayList<Double> aggregatedspeeds = new ArrayList<Double>();

        for(int j =0; j<coordinates.length;j++){
            int[] point = coordinates[j];
            int x = getXcoord(point[0]);
            int y = getYcoord(point[1]);
            /**
             * add to aggregatedpoints if not there
             */
            boolean already = false;
            int loc = -1;
            while(!already && loc<aggregatedpoints.size()-1){
                // check if the point exists already
                loc++;
                int[] xy = aggregatedpoints.get(loc);
                if(xy[0]==x && xy[1] == y){
                    already = true;
                }

            }

            if(!already){
                aggregatedpoints.add(new int[]{x,y});
                aggregatedspeeds.add(speeds[j]);
                counts.add(1);
            }else{
                double count = (double) counts.get(loc);
                aggregatedspeeds.set(loc, speeds[j]/(count+1)+aggregatedspeeds.get(loc)*count/(count+1));
                counts.set(loc,counts.get(loc)+1 );

            }
        }
        for(int j=0;j<aggregatedpoints.size();j++){
            int x = aggregatedpoints.get(j)[0];
            int y = aggregatedpoints.get(j)[1];
            double val = aggregatedspeeds.get(j);
            krigingLattice.addObservation(x,y,val);
        }
    }


    /**
     *  Current mean field
     * @return double 2d array
     */
    public double[][] getCurrentMean(){
        return this.currentMean;
    }

    /**
     * Infer the speed everywhere by running predict()
     * @return speed array using interpolation
     */
    public float[][] inferFullField(){
        int width1 = (int) Math.ceil(1.0*width/xreduction);
        int height1 = (int) Math.ceil(1.0*height/yreduction);
        double[][] coarsefield = new double[width1][height1];
        float[][] flcoarsefield = new float[width1][height1];
        boolean[][] decided = new boolean[width1][height1];

        for(int y=0;y<height1;y++){
            for(int x=0;x<width1;x++){
                coarsefield[x][y] = priorMeanSpeed;
                decided[x][y] = false;
            }
        }
        int bandwidth = 40;
        ArrayList<int[]> datapoints = krigingLattice.getPoints();
        boolean[] mask = new boolean[datapoints.size()];


        for(int position = 0; position<wavePositions.size();position++){
            for(int j=0; j<mask.length;j++){
                int[] xy = datapoints.get(j);
                if(Math.abs(wavePositions.get(position).get(invertXcoord(xy[0]),invertYcoord(xy[1])))<bandwidth){
                    mask[j] = true;
                }else mask[j] = false;
            }
            ArrayList<double[][]> inferred = krigingLattice.infer(mask);
            ArrayList<int[]> newpoints = new ArrayList<int[]>();
            for(int y=0;y<height1;y++){
                for(int x=0;x<width1;x++){
                    if(!decided[x][y]){
                        // check that we want to use this point
                        if(wavePositions.get(position).in(invertXcoord(x),invertYcoord(y))){
                            newpoints.add(new int[]{x,y});
                        }
                    }
                }
            }
            if(newpoints.size()>0){
                ArrayList<double[][]> predicted = krigingLattice.predict(inferred,newpoints);
                double[][] predictions = predicted.get(0);
                for(int j=0; j<predictions.length;j++){
                    int[] xy = newpoints.get(j);
                    if(!decided[xy[0]][xy[1]]) {
                        coarsefield[xy[0]][xy[1]] = predictions[j][0];
                        flcoarsefield[xy[0]][xy[1]] = (float) predictions[j][0];
                        decided[xy[0]][xy[1]] = true;

                    }
                }
            }

        }

        ImageProcessor flcoarse = new FloatProcessor(flcoarsefield);
        ImagePlus flimage = new ImagePlus("",flcoarse);
        flimage.show();
        flimage.updateAndDraw();


        double[] xvals = new double[width1];
        double[] yvals = new double[height1];
        for(int x=0;x<width1;x++){xvals[x] = x*xreduction;}
        for(int y=0;y<height1;y++){yvals[y] = y*yreduction;}
        int J = wavePositions.size();


        PiecewiseBicubicSplineInterpolator interpolator = new PiecewiseBicubicSplineInterpolator();
        PiecewiseBicubicSplineInterpolatingFunction	 poly = interpolator.interpolate(xvals,yvals,coarsefield);

        float[][] fullfield = new float[width][height];
        for(int y=0; y<height; y++){
            for(int x=0;x<width; x++){
                try {
                    fullfield[x][y] = (float) poly.value(x, y);
                }catch (OutOfRangeException e){
                    fullfield[x][y] = (float) priorMeanSpeed;
                }
            }
        }
        float[][] kernel = gaussianKernelWindow(1,3,12);
        FourierConvolution fc = new FourierConvolution();
        float[][] convolved = fc.convolvefft(fullfield,kernel);

        return convolved;
    }

    public double[][] gradNorm(double[][] function){
        double[][] newfunction = new double[function.length][function[0].length];

        return newfunction;
    }


    public float[][] listArrivalTimes(){
        ArrayList<Integer> xcoords = new ArrayList<Integer>();
        ArrayList<Integer> ycoords = new ArrayList<Integer>();
        ArrayList<Double> times = new ArrayList<Double>();
        for(int j=1; j<wavePositions.size();j++){
            int[][] boundaryCoordinates = wavePositions.get(j).getBoundaryCoordinates();
            for(int k=0; k<boundaryCoordinates.length;k++){
                int y = boundaryCoordinates[k][1];
                int x = boundaryCoordinates[k][0];
                if(x<=1||y<=1 || x>= width-2 || y>= height-2) continue;
                ycoords.add(y);
                xcoords.add(x);
                times.add((double)this.times.get(j));
            }
        }
        float[][] out = new float[xcoords.size()][3];
        for(int i=0;i<xcoords.size();i++){
            out[i][0] = xcoords.get(i);
            out[i][1] = ycoords.get(i);
            out[i][2] = times.get(i).floatValue();
        }
        return out;

    }

    public float[][] iterateSpeedInference(float[][] initialguess,float[][] kernel,int iterations){
        FourierConvolution fc = new FourierConvolution();
        float[][] lambda = fc.convolvefft(initialguess,kernel);
        float[][] C = new float[width][height];
        float[][] T = new float[width][height];

        for(int y=0;y<height;y++){
            for(int x=0; x<width;x++){
                C[x][y] = initialguess[x][y];
            }
        }
        for(int i=0;i<iterations;i++){
            lambda = fc.convolvefft(C,kernel);
        }
        return null;
    }


    /**
     * Generate a spherical gaussian window of width with height and spread parameters
     * @param height
     * @param spread
     * @param width
     * @return
     */
    public float[][] gaussianKernelWindow(float height, float spread, int width){
        float[][] kernel = new float[width*2+1][width*2+1];
        double prefactor = 1.0/Math.sqrt(2*Math.PI)/spread*height;
        for(int y=0;y<2*width+1;y++){
            for(int x=0;x<2*width+1;x++){
                kernel[x][y] = (float)( prefactor*Math.exp(-Math.sqrt((x-width)*(x-width) + (y-width)*(y-width))/2/spread));
            }
        }
        return kernel;
    }

    public ArrayList<double[][]> sampleSpeedFields(int numSamples, ImplicitShape2D shape, double bandwidth){
        ArrayList<int[]> datapoints = krigingLattice.getPoints();

        boolean[] mask = new boolean[datapoints.size()];

        /**
         * Figure out which lattice points are within our desired region within shape +- bandwidth
         */

        for(int i=0;i<mask.length;i++){
            int[] xy = datapoints.get(i);
            if(Math.abs(shape.get(invertXcoord(xy[0]),invertYcoord(xy[1])))<bandwidth){
                mask[i] = true;
            }else mask[i] = false;
        }
        ArrayList<double[][]> inferred = krigingLattice.infer(mask);

        ArrayList<int[]> newpoints = new ArrayList<int[]>();
        ArrayList<Integer> whichpoint = new ArrayList<Integer>(); // which point in newpoints does this x,y location correspond to?

        int num_new_points = 0;
        for(int y=0;y<height;y++){
            for(int x=0;x<width;x++){
                if(shape.get(x,y)>-bandwidth && shape.get(x,y)<0){
                    // check first to make sure that this point is not already in newpoints
                    int x2 = getXcoord(x);
                    int y2 = getYcoord(y);
                    /**
                     * add to newpoints if not there
                     */
                    boolean already = false;
                    int loc = -1;
                    while(!already && loc<newpoints.size()-1){
                        // check if the point exists already
                        loc++;
                        int[] xy = newpoints.get(loc);
                        if(xy[0]==x2 && xy[1] == y2){
                            already = true;
                        }

                    }
                    if(already) {
                        // point is already in
                        whichpoint.add(loc);

                    }else{
                        whichpoint.add(num_new_points);
                        newpoints.add(new int[]{x2,y2});
                        num_new_points++;
                    }

                }
            }
        }

        /**
         * Figure out which of the data points we will use
         */

        if(datapoints.size()==0 || newpoints.size()==0){
            ArrayList<double[][]> out = new ArrayList<double[][]>(numSamples);
            // sample from the prior
            // priormean + LIZ  (Z is random vector)
            for(int i=0;i<numSamples;i++){
                // solve for convolution
                double[][] speed = new double[width][height];
                for(int y=0;y<height;y++){
                    for(int x=0; x<width; x++){
                        speed[x][y] = Math.abs(this.priorMeanSpeed + Normal.staticNextDouble(0,Math.sqrt(this.priorVarSpeed)));
                    }
                }
                out.add(speed);
            }
            return out;
        }

        /**
         * Pull samples out of this.krigingLattice, and convert those samples to speed fields
         */


        ArrayList<double[][]> predicted = krigingLattice.predict(inferred,newpoints);
        ArrayList<double[][]> krigingsamples = krigingLattice.sample(predicted,numSamples);
        ArrayList<double[][]> samples = new ArrayList<double[][]>(numSamples+8);

        for(int i=0; i<numSamples; i++){

            samples.add(new double[width][height]);
            int ptcount = 0;
            for(int y=0; y<height; y++){
                for(int x=0; x<width; x++){
                    if(shape.get(x,y)>-bandwidth && shape.get(x,y)<0){
                        samples.get(i)[x][y] = krigingsamples.get(i)[whichpoint.get(ptcount)][0];
                        ptcount++;
                    }
                    else{
                        samples.get(i)[x][y] = inferred.get(0)[0][0]+ Normal.staticNextDouble(0,Math.sqrt(inferred.get(1)[0][0]));
                    }
                }
            }
        }
        for(int i=0;i<8;i++){
            double[][] speed = new double[width][height];
            for(int y=0;y<height;y++){
                for(int x=0; x<width; x++){
                    speed[x][y] = Math.abs(this.priorMeanSpeed + Normal.staticNextDouble(0,Math.sqrt(this.priorVarSpeed)));
                }
            }
            samples.add(speed);
        }


        return samples;
    }



    public int getXcoord(int x){
        return (int) Math.round(x/this.xreduction);
    }

    public int invertXcoord(int x){ return (int) x*this.xreduction; }

    public int getYcoord(int y){
        return (int)Math.round(y/this.yreduction);
    }

    public int invertYcoord(int y){ return (int) y*this.yreduction;}
}
