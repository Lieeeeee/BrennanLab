package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class Cumulative_sum implements PlugInFilter {
    ImagePlus imp;

    @Override
    public void run(ImageProcessor ip) {
        int slices = imp.getStackSize();
        int width = imp.getWidth();
        int height = imp.getHeight();
        ImageStack newStack = imp.createEmptyStack();
        float[] pixels;
        float[] priorpixels = new float[width * height];
        for (int j = 0; j < priorpixels.length; j++) {
            priorpixels[j] = 0;
        }
        FloatProcessor newProcessor;
        for (int i = 1; i <= slices; i++) {
            IJ.showProgress(i, slices);
            pixels = new float[width * height];
            pixels = (float[]) imp.getStack().getProcessor(
                    i <= slices ? i : slices).getPixels();
            for (int j = 0; j < priorpixels.length; j++) {
                priorpixels[j] += pixels[j];
            }
            newProcessor = new FloatProcessor(width, height, priorpixels
                    .clone(), ip.getColorModel());
            newStack.addSlice(imp.getStack().getSliceLabel(i), newProcessor);

        }
        ImagePlus newImage = new ImagePlus(imp.getTitle(), newStack);
        newImage.show();
        newImage.updateAndDraw();
        newImage.unlock();
        imp.updateAndDraw();
        imp.unlock();
    }

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return STACK_REQUIRED + DOES_32;
    }

}