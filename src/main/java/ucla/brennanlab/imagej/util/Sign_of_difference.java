package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class Sign_of_difference implements PlugInFilter {
    ImagePlus imp; // input image
    ImagePlus difference;
    String arg;
    int N;

    public float mean(float[] nums) {
        float sum = 0;
        for (int j = 0; j < nums.length; j++) {
            sum += nums[j];

        }
        return (sum / nums.length);
    }

    public void run(ImageProcessor ip) {
        byte[] pixels;
        byte[] prevpixels;
        float[] differences;
        int width = ip.getWidth();
        int height = ip.getHeight();
        ImageStack stack = imp.getStack();
        ImageStack sgnstack = imp.createEmptyStack();
        prevpixels = new byte[width * height];
        prevpixels = (byte[]) stack.getPixels(1);

        FloatProcessor differenceIP;
        for (int i = 1; i < N; i++) {
            differences = new float[width * height];
            pixels = (byte[]) stack.getPixels(i + 1);
            for (int j = 0; j < width * height; j++) {
                differences[j] = sgn(pixels[j] - prevpixels[j]);
                prevpixels[j] = pixels[j];
                IJ.showProgress(i, N);
            }
            differenceIP = new FloatProcessor(width, height, differences, null);
            sgnstack.addSlice("" + mean(differences), differenceIP);
        }
        difference = new ImagePlus("Sign of Differences", sgnstack);
        difference.updateAndDraw();

        difference.show();
    }

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.N = imp.getStackSize();
        return STACK_REQUIRED + NO_CHANGES + DOES_8G;

    }

    public byte sgn(int num) {
        if (num > 0)
            return 1;
        else if (num < 0)
            return -1;
        else
            return 0;

    }
}