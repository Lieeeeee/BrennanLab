package ucla.brennanlab.imagej.util.levelsets.tests;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;
import ucla.brennanlab.imagej.util.levelsets.fastMarchingSolver;

public class levelsetTest implements PlugIn {

    @Override
    public void run(String arg) {
        int width = 320;
        int height = 240;
        boolean[][] mask = new boolean[width][height];
        float[][] speeds = new float[width][height];
        RandomEngine engine = new cern.jet.random.engine.MersenneTwister(
                new java.util.Date());
        Normal normalRV1 = new Normal(4, 6, engine);
        Normal normalRV2 = new Normal(5, 6, engine);

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                mask[i][j] = Math.sqrt(Math.pow(i - width / 2, 2)
                        + Math.pow(j - height / 2, 2)) < 50 ? true : false;
                speeds[i][j] = (float) ((float) j > height / 2 ? Math.max(normalRV1
                        .nextDouble(), 1) : Math.max(normalRV2.nextDouble(), 1)); // test
                // unit
                // speed
            }
        }
        ImplicitShape2D ls = new ImplicitShape2D(mask);
        Roi r = ls.getRoi(false);
        /*
		ImagePlus imp = ls.getLSImagePlus();
		imp.setRoi(r);
		imp.setTitle("levelSetSignedDistance");
		imp.show();
		imp.updateAndDraw();
		imp.getProcessor().setRoi(r);
		/**
		 * Test fast marching
		 */

        IJ.log("Testing fast marching propagation...");
        long startTime = System.nanoTime();
        long endTime;
        fastMarchingSolver fm;
        ImplicitShape2D ri;
        try {
            fm = new fastMarchingSolver(ls, speeds);
            ri = fm.solveAndReturnImplicitShape2D(20);
        } finally {
            endTime = System.nanoTime();

        }
        IJ.log("Finished fast-marching in " + (float) (endTime - startTime)
                / 1000000000 + " s");
        ImagePlus imp1 = ri.getLSImagePlus();
        imp1.setTitle("propagated");
        imp1.setRoi(ri.getRoi(false));
        Overlay ov = new Overlay(ri.getRoi(false));
        imp1.setOverlay(ov);
        imp1.show();
        imp1.updateAndDraw();
        imp1.getProcessor().draw(r);

    }

}