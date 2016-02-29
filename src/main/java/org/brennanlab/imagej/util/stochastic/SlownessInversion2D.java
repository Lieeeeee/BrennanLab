package org.brennanlab.imagej.util.stochastic;
import org.brennanlab.imagej.util.FourierConvolution;

import java.util.ArrayList;

/**
 * Created by changjc on 2/28/16.
 */
public class SlownessInversion2D {
    private int width, height;
    public ArrayList<int[]> points;
    public ArrayList<Double> times;
    private float[][] R; // regularization kernel
    private float[][] observationField;
    private double sigma_spatial = 2.0;
    private double sigma_temporal = 0.25;
    private float[][] partialYKernel = new float[][]{ {0,-1,0}, {0,2,0}, {0,-1,0}};
    private float[][] partialXKernel = new float[][]{ {0,0,0},{-1,2,-1},{0,0,0}};
    private float[][] laplacianKernel = new float[][] {};

    private float[][] observationKernel;
    FourierConvolution fft;

    SlownessInversion2D(int width, int height){
        this.width = width;
        this.height = height;
        FourierConvolution fft = new FourierConvolution();
        this.observationField = new float[width][height];
        int kw = (int) Math.ceil(this.sigma_spatial * 3);
        observationKernel = new float[kw*2+1][kw*2+1];
        for(int y=-kw;y<=kw;y++)
            for(int x=-kw;x<=kw;x++){
                observationKernel[x+kw][y+kw] = (float) (Math.exp( -(x-y)*(x-y)/2/sigma_spatial)/Math.sqrt(2*Math.PI*sigma_spatial)/sigma_temporal/sigma_temporal );
            }

        return;
    }
    public void addArrivalTimes(ArrayList<int[]> points, ArrayList<Double> times){
        if(times.size() != points.size()) throw new IllegalArgumentException("There should be  an arrival time for each point added");
        this.points.addAll(points);
        this.times.addAll(times);
    }

    public float[][] generateObservationField(ArrayList <int[]> points, ArrayList<Double> times, float[][] T){
        // For each ArrivalTime, add Gaussian
        float[][] observationField = new float[width][height];
        int kw = (observationKernel.length-1)/2;
        for(int j=0; j<points.size();j++){
            int x0 = points.get(j)[0];
            int y0 = points.get(j)[1];
            // Explore the window of values for T around this point, and multiply the residuals by the observation kernel
            for(int y = y0-kw; y<= y0+kw; y++){
                for(int x = x0-kw; x<- x0+kw; x++){
                    if(x<0 || y<0 || x>=width || y>= height) continue;
                    else{
                        observationField[x][y] += observationKernel[x+kw-x0][y+kw-y0]*(times.get(j)-T[x][y]);
                    }
                }
            }
        }

        return observationField;
    }

    public ArrayList<float[][]> nextIterate(ArrayList<float[][]> currentState){
        if(currentState.size()!=3) throw new IllegalArgumentException("Need three elements for current state (\\lambda, T,C)");
        float[][] lambda = currentState.get(0);
        float[][] T = currentState.get(1);
        float[][] C = currentState.get(2);

        float[][] TNext = new float[width][height];
        float[][] CNext = new float[width][height];

        // Solve \lambda = RC
        float[][] lambdaNext = fft.convolvefft(C,R);

        // Solve Poisson equation using Born approximation

        float[][] g = generateObservationField(points,times,T);
        float[][] epsilon = new float[width][height];

        float[][] Tx = fft.convolvefft(T,partialXKernel);
        float[][] Ty = fft.convolvefft(T,partialXKernel);

        for(int y=0; y<height;y++){
            for(int x=0;x<width;x++){
                epsilon[x][y] = (float)(lambdaNext[x][y]/Math.sqrt(Tx[x][y]*Tx[x][y]+Ty[x][y]*Ty[x][y]) );
            }
        }

        float[][] source = new float[width][height];
        for(int y=0; y<height; y++){
            for(int x=0; x<width; x++){
                source[x][y] = g[x][y]/epsilon[x][y];
            }
        }

        // Use first born approximation to solve for T first
        TNext = fft.convolvefft(source,laplacianKernel);

        // Implement first-order correction for Born approximation

        // Solve the constraint C = |\nabla T|

        Tx = fft.convolvefft(TNext,partialXKernel);
        Ty = fft.convolvefft(TNext,partialXKernel);

        for(int y=0; y<height;y++){
            for(int x=0; x<width;x++){
                CNext[x][y] = (float) Math.sqrt(Tx[x][y]*Tx[x][y] + Ty[x][y]*Ty[x][y]);
            }
        }


        ArrayList<float[][]> nextState = new ArrayList<float[][]>(3);
        nextState.add(lambdaNext);
        nextState.add(TNext);
        nextState.add(CNext);
        return nextState;
    }
}
