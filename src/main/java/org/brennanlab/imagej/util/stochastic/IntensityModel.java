package org.brennanlab.imagej.util.stochastic;

import ij.ImageStack;
import ij.process.ImageProcessor;

public interface IntensityModel {
    void Infer(ImageProcessor ip,boolean[][] labels);
    void Infer(ImageProcessor ip, boolean[][] labels,boolean updatePrior);
    void Infer(ImageProcessor ip, boolean[][] labels,boolean updatePrior, int innerband);

    /**
     * For inference on a stack
     * @param is
     * @param labels
     */
    void Infer(ImageStack is, boolean[][][] labels);

    double logpOut(double pixelval);

    double logpIn(double pixelval);
    double getPosteriorMean(boolean in);
    double getPosteriorPrecision(boolean in);

    void Infer(ImageStack stack, boolean[][][] mask, boolean updatePrior);
    void multiplyPrecision(double d);

}