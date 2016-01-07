package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;

public class SegmentationError implements PlugInFilter {
    Roi segmentation;
    ImagePlus imp;

    public void run(ImageProcessor ip) {

        IJ
                .log("This is a tool to compare the segmentation given by an ROI against a ground truth");
        IJ
                .log("This tool calculates the mismatched pixels vs the true boundary length");

        boolean[][] mask = new boolean[ip.getWidth()][ip.getHeight()];
        for (int x = 0; x < ip.getWidth(); x++) {
            for (int y = 0; y < ip.getHeight(); y++) {
                mask[x][y] = ip.getPixelValue(x, y) == 0 ? true : false;
            }
        }

        ImplicitShape2D rls = new ImplicitShape2D(mask);

        int mismatched = 0;
        double len = 0;

        for (int x = 1; x < this.imp.getWidth() - 1; x++) {
            for (int y = 1; y < this.imp.getHeight() - 1; y++) {
                if (Math.abs(rls.get(x, y)) < 1) {
                    len += 1.0 / Math.PI / (rls.get(x, y) * rls.get(x, y) + 1);
                }
                if (mask[x][y] != this.segmentation.contains(x, y)) {
                    mismatched++;
                }

            }
        }
        double miss = 1.0 * mismatched / len;
        IJ.log("Mismatched " + miss);

        imp.unlock();

    }

    @Override
    public int setup(String arg, ImagePlus imp) {
        // TODO Auto-generated method stub
        this.segmentation = imp.getRoi();
        this.imp = imp;
        return ROI_REQUIRED + DOES_8G;
    }

}