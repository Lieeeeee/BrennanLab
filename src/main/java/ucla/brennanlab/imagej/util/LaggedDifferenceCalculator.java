package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class LaggedDifferenceCalculator implements PlugInFilter {
    ImagePlus imp; // input image
    ImagePlus difference;
    String arg;
    boolean filter;
    int N;
    int lag = 1;

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

        FloatProcessor differenceIP;
        for (int i = 1 + lag; i <= N; i++) {
            prevpixels = (byte[]) stack.getPixels(i - lag);
            differences = new float[width * height];
            pixels = (byte[]) stack.getPixels(i);
            for (int j = 0; j < width * height; j++) {
                differences[j] = pixels[j] - prevpixels[j];
                if (filter && (Math.abs(differences[j]) > 20))//TODO use Huber function or something
                    differences[j] = 0;
            }

            differenceIP = new FloatProcessor(width, height, differences, null);

            sgnstack.addSlice("" + mean(differences), differenceIP);
            IJ.showProgress(i, N);
        }

        difference = new ImagePlus("Differences", sgnstack);
        difference.updateAndDraw();

        difference.show();
    }

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.N = imp.getStackSize();
        this.filter = false;
        if (!(arg.length() > 0)) {
            arg = Macro.getOptions();
        }

        if (arg != null && !arg.isEmpty()) {
            String[] args = arg.split(";");
            for (int i = 0; i < args.length; i++) {
                String[] param = args[i].split("=");
                if (param[0].equals("filter")) {
                    //IJ.showMessage("filtering");
                    filter = true;
                } else if (param[0].equals("lag")) {
                    lag = Integer.parseInt(param[1].trim());
                    //IJ.showMessage("Using lag of "+lag);
                }
            }
        }
        return STACK_REQUIRED + NO_CHANGES + DOES_8G + DOES_16 + DOES_32;

    }
}