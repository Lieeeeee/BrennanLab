package org.brennanlab.imagej.csd;

import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class SelectionToText implements PlugInFilter {
    ImagePlus imp;
    Roi r;

    protected GenericDialog params() {
        GenericDialog gd = new GenericDialog("Parameters");
        gd.addNumericField("Width divisions", 10, 0);
        gd.addNumericField("Height divisions", 8, 0);
        gd.addCheckbox("Image", true);
        return gd;
    }

    @Override
    public void run(ImageProcessor ip) {
        GenericDialog gd = params();
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        int widthdiv = (int) gd.getNextNumber();
        int heightdiv = (int) gd.getNextNumber();
        boolean img = gd.getNextBoolean();
        int[][] pixeltotal = new int[widthdiv][heightdiv];
        int[][] incounts = new int[widthdiv][heightdiv];

        for (int y = 0; y < heightdiv; y++) {
            for (int x = 0; x < widthdiv; x++) {
                pixeltotal[x][y] = 0;
                incounts[x][y] = 0;
            }
        }

        for (int y = 0; y < ip.getHeight(); y++) {
            int py = (int) (Math.floor(1.0 * y * heightdiv / ip.getHeight()));
            for (int x = 0; x < ip.getWidth(); x++) {

                int px = (int) (Math.floor(1.0 * x * widthdiv / ip.getWidth()));
                pixeltotal[px][py]++;
                if (r.contains(x, y))
                    incounts[px][py]++;
            }
        }

        ResultsTable rt = new ResultsTable();

        for (int y = 0; y < heightdiv; y++) {
            rt.incrementCounter();
            for (int x = 0; x < widthdiv; x++) {

                rt.addValue("" + x, (float) incounts[x][y] / pixeltotal[x][y]);
            }
        }

        if (img) {
            FloatProcessor ip1 = new FloatProcessor(ip.getWidth(), ip
                    .getHeight());
            for (int x = 0; x < ip.getWidth(); x++) {
                for (int y = 0; y < ip.getHeight(); y++) {
                    ip1.setf(x, y, (float) incounts[(int) (Math.floor(1.0 * x
                            * widthdiv / ip.getWidth()))][(int) (Math.floor(1.0
                            * y * heightdiv / ip.getHeight()))]
                            / pixeltotal[(int) (Math.floor(1.0 * x * widthdiv
                            / ip.getWidth()))][(int) (Math.floor(1.0
                            * y * heightdiv / ip.getHeight()))]);
                }
            }
            ImagePlus im = new ImagePlus("", ip1);
            im.show();
            im.updateAndDraw();
        }

        rt.show("Binned ROI region occupations");

    }

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.r = imp.getRoi();
        return ROI_REQUIRED + DOES_ALL;
    }

}