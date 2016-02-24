package org.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.Arrays;

public class Change_from_baseline implements PlugInFilter {
    public static final String[] baseMethods = {"Mean", "Median", "Max", "Min"};
    public static final int MEAN_METHOD = 0;
    public static final int MEDIAN_METHOD = 1;
    public static final int MAX_METHOD = 2;
    public static final int MIN_METHOD = 3;
    public static final String[] scaleMethods = {"None", "Baseline mean",
            "Baseline Standard Deviation", "Difference in Absorbance"};
    public static final int NO_SCALE = 0;
    public static final int MEAN_SCALE = 1;
    public static final int SD_SCALE = 2;
    public static final int LOG_SCALE = 3;
    public static final String[] PROJECTIONS = {"Mean", "Median", "SD", "Max",
            "Min"};
    public static final int MEAN = 0;
    public static final int MEDIAN = 1;
    public static final int SD = 2;
    public static final int MAX = 3;
    public static final int MIN = 4;
    public static int basemethod = MEAN_METHOD;
    public static int scalemethod = NO_SCALE;
    ImagePlus original;
    ImagePlus baseline;
    ImagePlus baselineProjection;
    ImagePlus resultBlob;
    String arg;
    int baseStart = 2, baseStop = 9;

    protected GenericDialog OptionSelector() {
        GenericDialog gd = new GenericDialog("Options", IJ.getInstance());
        gd.addMessage("Select from the options below:");
        gd.addNumericField("Baseline start slice:", baseStart, 0);
        gd.addNumericField("stop slice:", baseStop, 0);
        gd.addChoice("Projection", baseMethods, baseMethods[basemethod]);
        gd.addChoice("Result Scaling method", scaleMethods,
                scaleMethods[scalemethod]);
        return (gd);

    }

    public void run(ImageProcessor ip) {
        GenericDialog gd = OptionSelector();
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        baseStart = (int) gd.getNextNumber();
        baseStop = (int) gd.getNextNumber();
        basemethod = gd.getNextChoiceIndex();
        scalemethod = gd.getNextChoiceIndex();

        original = WindowManager.getCurrentImage();

        if (original == null) {
            IJ.noImage();
            return;
        }

        ImageStack stack = original.getStack();
        ImageStack baselineStack = original.createEmptyStack();

        for (int i = baseStart; i <= baseStop; i++) {

            ImageProcessor ip1 = stack.getProcessor(i <= baseStop ? i
                    : baseStop);
            ip1 = ip1.crop();
            baselineStack.addSlice("", ip1);

        }
        baseline = new ImagePlus("", baselineStack);

        baselineProjection = zProjector(baseline, basemethod);

        // baselineProjection.show();
        // baselineProjection.updateAndDraw();

        subtractAndScale();
        resultBlob.show();
        resultBlob.updateAndDraw();

    }

    public int setup(String arg, ImagePlus img) {
        this.arg = arg;
        this.original = img;
        return STACK_REQUIRED + DOES_8G + DOES_8C + DOES_16 + DOES_32
                + NO_CHANGES + NO_UNDO;
    }

    private void subtractAndScale() {
        // Just use ImagePlus' subtract function
        ImageStack resultstack = original.createEmptyStack();

        int numSlices = original.getStackSize();
        double difference;

        ImageProcessor ip1;
        ImageProcessor ipBase = baselineProjection.getProcessor();
        FloatProcessor ip32;
        int width = original.getWidth();
        int height = original.getHeight();
        int dimension = width * height;
        float[] basepixels = new float[dimension];
        basepixels = (float[]) ipBase.getPixels();

        float[] scalepixels = new float[dimension];
        scalepixels = basepixels;

        float[][] pixels = new float[width][height];
        ImagePlus scalebase = new ImagePlus();
        if (scalemethod != LOG_SCALE) {
            scalebase = zProjector(baselineProjection, scalemethod);
            scalepixels = (float[]) scalebase.getProcessor().getPixels();
        }

        for (int i = 1; i <= numSlices; i++) {
            // ip1=new ImageProcessor();
            IJ.showProgress(i, numSlices);
            ip1 = original.getStack().getProcessor(
                    i <= numSlices ? i : numSlices);
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    if (scalemethod == LOG_SCALE) {
                        difference = Math.log(ip1.getPixelValue(x, y) + 1)
                                - Math.log(basepixels[y * width + x] + 1);
                        difference = difference / Math.log(10);
                    } else {
                        difference = ip1.getPixelValue(x, y)
                                - basepixels[y * width + x];

                        if (scalemethod == MEAN_SCALE)
                            difference /= scalepixels[y * width + x];
                    }
                    pixels[x][y] = (float) difference;
                }
            }
            ip32 = new FloatProcessor(pixels);
            resultstack.addSlice("" + i, ip32);
        }

        resultBlob = new ImagePlus("Scaled by " + scaleMethods[scalemethod],
                resultstack);
        return;

    }

    public ImagePlus zProjector(ImagePlus img, int METHOD) {
        // first get pixels into 3d pixel stack array
        int nSlices = img.getStackSize();
        if (nSlices < 2)
            return img;
        int height = img.getHeight();
        int width = img.getWidth();
        float[][][] pixels = new float[nSlices][width][height];
        ImageProcessor ip;
        for (int slice = 1; slice <= nSlices; slice++) {
            ip = img.getStack().getProcessor(slice);
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    pixels[slice - 1][x][y] = ip.getPixelValue(x, y);
                }
            }

        }
        float[][] projectedpixels = new float[width][height];
        float[] zPixels = new float[nSlices];
        float resultPixel = 0;
        double sum1, sum2, sd = 0;
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int slice = 0; slice < nSlices; slice++) {
                    zPixels[slice] = pixels[slice][x][y];
                }
                switch (METHOD) {
                    case SD:
                        sum1 = 0;
                        sum2 = 0;
                        for (int i = 0; i < zPixels.length; i++) {
                            sum1 += zPixels[i];
                            sum2 += (Math.pow(zPixels[i], 2));
                        }
                        sd = Math.sqrt((nSlices * sum2 - Math.pow(sum1, 2))
                                / (nSlices * (nSlices - 1)));
                        // IJ.showMessage(""+sd);
                        resultPixel = (float) sd;
                        break;
                    case MEAN:
                        sum1 = 0;
                        for (int i = 0; i < zPixels.length; i++) {
                            sum1 += zPixels[i];
                        }
                        resultPixel = (float) (sum1 / nSlices);
                        break;
                    case MEDIAN:
                        Arrays.sort(zPixels);
                        if (nSlices % 2 == 0)
                            resultPixel = (zPixels[nSlices / 2] + zPixels[nSlices / 2 + 1]) / 2;
                        else
                            resultPixel = zPixels[(int) (Math
                                    .floor(nSlices / 2.0) + 1)];
                        resultPixel = 0;
                    case MIN:
                        Arrays.sort(zPixels);
                        resultPixel = zPixels[0];
                        break;
                    case MAX:
                        Arrays.sort(zPixels);
                        resultPixel = zPixels[nSlices - 1];
                        break;
                }
                projectedpixels[x][y] = resultPixel;
            }
        }
        FloatProcessor ip2 = new FloatProcessor(projectedpixels);
        ImagePlus projected = new ImagePlus("zProjection", ip2);
        return projected;
    }

}
