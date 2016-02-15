package ucla.brennanlab.imagej.util.stochastic;

import cern.jet.random.Normal;
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
    public double betahat; // Temporary
    public double Abeta; // Temporary
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
        this.betahat = priormean;
        this.standardNormal = new Normal(0, 1,
                new cern.jet.random.engine.MersenneTwister(new java.util.Date())
        );
        krigingLattice = new Kriging2DLattice(priormean,1.0/priorvar,width/xreduction,height/yreduction);
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

        addObservations(incomingBoundaryCoordinates,incomingBoundarySpeeds);
        this.times.add(position,time);
        this.wavePositions.add(position,incomingShape);

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
            while(!already && loc<aggregatedpoints.size()){
                // check if the point exists already
                loc++;
                int[] xy = aggregatedpoints.get(loc);
                if(xy[0]==x && xy[1] == y){
                    already = true;
                }

            }

            if(!already){
                aggregatedpoints.add(point);
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
     *
     * @return
     */
    public double[][] getCurrentMean(){
        return this.currentMean;
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

    public ArrayList<double[][]> sampleSpeedFields(int numSamples, ImplicitShape2D shape, double bandwidth){
        ArrayList<int[]> datapoints = krigingLattice.getPoints();

        if(datapoints.size()==0){
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

        boolean[] mask = new boolean[datapoints.size()];

        /**
         * Figure out which lattice points are within our desired region within shape +- bandwidth
         */

        for(int i=0;i<mask.length;i++){
            int[] xy = datapoints.get(i);
            if(Math.sqrt(Math.abs(shape.get(invertXcoord(xy[0]),invertYcoord(xy[1]))))<bandwidth){
                mask[i] = true;
            }else mask[i] = false;
        }
        ArrayList<double[][]> inferred = krigingLattice.infer(mask);

        ArrayList<int[]> newpoints = new ArrayList<int[]>();

        /**
         * Figure out which of the data points we will use
         */



        /**
         * Pull samples out of this.krigingLattice, and convert those samples to speed fields
         */


        ArrayList<double[][]> predicted = krigingLattice.predict(inferred,newpoints);
        ArrayList<double[][]> krigingsamples = krigingLattice.sample(predicted,numSamples);

        ArrayList<double[][]> samples = new ArrayList<double[][]>(numSamples);

        for(int i=0; i<numSamples; i++){

            samples.add(new double[width][height]);
            for(int x=0; x<width; x++){
                for(int y=0; y<height; y++){
                    if(shape.get(x,y)>bandwidth*bandwidth){
                        samples.get(i)[x][y] = this.betahat;
                    }
                    else{
                        // interpolate

                    }
                }
            }
        }


        return samples;
    }



    public int getXcoord(int x){
        return (int) Math.floor(x/this.xreduction);
    }

    public int invertXcoord(int x){ return (int) x*this.xreduction; }

    public int getYcoord(int y){
        return (int)Math.floor(y/this.yreduction);
    }

    public int invertYcoord(int y){ return (int) y*this.yreduction;}
}
