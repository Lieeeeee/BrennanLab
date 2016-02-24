package org.brennanlab.imagej.util;

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import org.jtransforms.fft.FloatFFT_2D;

import java.util.ArrayList;

/**
 * Created by changjc on 2/24/16.
 */
public class FourierConvolution implements PlugInFilter {
    ImagePlus imp;
    String arg;
    public void run(ImageProcessor ip) {


    }

    public FourierConvolution(){

    }

    public ArrayList<float[][]> pad(float[][] original){
        //int paddimensionx = (int) ;
        return null;
    }

    public float[][] dePad(float[][] padded,int xstart, int ystart, int width, int height){
        float[][] unpadded = new float[width][height];
        for(int y=0;y<height;y++){
            for(int x=0;x<width;x++){
                unpadded[x][y] = padded[x+xstart][y+ystart];
            }
        }

        return unpadded;
    }

    public float[][] convolvefft(float[][] function, float[][] kernel){
        int paddimension = (int) Math.pow(2,
                (int) Math.ceil(Math.log(3*Math.max(Math.max(function.length,function[0].length)
                        ,Math.max(kernel.length,kernel[0].length) ))/Math.log(2) ));
        float[][] paddedfunction = new float[paddimension][paddimension];
        float[][] paddedkernel = new float[paddimension][paddimension];

        int xshift = paddimension/2-function.length/2;
        int yshift = paddimension/2-function[0].length/2;


        for(int y=0;y<function[0].length;y++){
            for(int x=0;x<function.length;x++){
                paddedfunction[x+xshift][(y+yshift)] = function[x][y];
            }
        }


        int kernelxshift = paddimension/2-kernel.length/2;
        int kernelyshift = paddimension/2-kernel[0].length/2;

        for(int y=0;y<kernel[0].length;y++){
            for(int x=0;x<kernel.length;x++){
                paddedkernel[x+kernelxshift][(y+kernelyshift)] = kernel[x][y];
            }
        }


        /**
         * Pad the function and the kernel
         */
        FloatFFT_2D fft = new FloatFFT_2D(paddimension,paddimension);
        fft.realForward(paddedkernel);
        fft.realForward(paddedfunction);

        float[][] paddedconvolved = new float[paddimension][paddimension];
        for(int y=0;y<paddimension;y++){
            for(int x=0;x<paddimension;x++){
                paddedconvolved[x][y] = paddedkernel[x][y]*paddedfunction[x][y];
            }
        }

        fft.realInverse(paddedconvolved,true);

        float[][] convolved = new float[function.length][function[0].length];
        for(int y=0;y<function[0].length;y++){
            for(int x=0;x<function.length;x++){
                convolved[x][y] = paddedconvolved[x+xshift][(y+yshift)];
            }
        }
        return paddedconvolved;
    }


    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.arg = arg;
        return DOES_ALL + NO_CHANGES;
    }
}
