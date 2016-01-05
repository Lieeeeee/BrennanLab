package ucla.brennanlab.imagej.util.levelsets.tests;

/**
 * This class advects an image from a user-specified ROI or point 
 * according to speed given as image, at 1Hz, and generates a CSD-like
 * image sequence
 */
import ucla.brennanlab.imagej.util.levelsets.GenericLevelSet;
import ucla.brennanlab.imagej.util.levelsets.RoiImageLevelSet;
import ucla.brennanlab.imagej.util.levelsets.fastMarchingSolver;
import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class fastMarchingTest implements PlugInFilter {

    ImagePlus imp;
    final float imageNoise = (float) 4; // sd of image noise
    final float outerMean = 0;
    final float innerMean = 3;
    final int time = 300; // number of frames to do;

    RandomEngine engine = new cern.jet.random.engine.MersenneTwister(
	    new java.util.Date());
    Normal normalRV1 = new Normal(0, 0.2, engine); // speed noise
    Normal normalInner = new Normal(innerMean, imageNoise, engine);
    Normal normalOuter = new Normal(outerMean, imageNoise, engine);

    public void run(ImageProcessor ip) {
	Roi r = imp.getRoi();

	int width = ip.getWidth();
	int height = ip.getHeight();
	float[][] speeds = new float[width][height];

	for (int i = 0; i < ip.getWidth(); i++) {
	    for (int j = 0; j < ip.getHeight(); j++) {
		speeds[i][j] = (float) Math.max(ip.getPixelValue(i, j), 0.1);
		// + normalRV1.nextDouble(), 1);
	    }
	}
	int t = 1;
	ImagePlus result = imp.createImagePlus();
	ImagePlus groundtruth = imp.createImagePlus();
	groundtruth.setTitle("Ground_truth");
	result.setTitle("Synthetic_image_stack");
	IJ.log("Testing fast marching propagation...");
	long startTime = System.nanoTime();
	long endTime;
	fastMarchingSolver fm = new fastMarchingSolver(new GenericLevelSet(
		RoiImageLevelSet.roiToMask(r, width, height)), speeds);
	ImageStack syntheticWaveStack = imp.createEmptyStack();
	ImageStack syntheticGroundTruthStack = imp.createEmptyStack();
	syntheticWaveStack.addSlice("UCLA " + t, this
		.maskToProcessor(RoiImageLevelSet.roiToMask(r, width, height)));
	    syntheticGroundTruthStack.addSlice("Ground_Truth_" + t, this
		    .maskToByteProcessor(RoiImageLevelSet.roiToMask(r, width, height)));	
	GenericLevelSet gls;

	while (t < time) {
	    IJ.showProgress(t, time);
	    IJ.log("Generating image at time " + t);
	    try {
		gls = fm.solveAndReturnSignedDistance(t);
	    } finally {
		endTime = System.nanoTime();

	    }
	    IJ.log("Finished fast-marching in " + (float) (endTime - startTime)
		    / 1000000000 + " s");
	    syntheticWaveStack.addSlice("Synthetic_" + t, this
		    .maskToProcessor(gls.getMask()));
	    syntheticGroundTruthStack.addSlice("Ground_Truth_" + t, this
		    .maskToByteProcessor(gls.getMask()));

	    if (fm.isDone())
		break;
	    t++;
	    result.setStack(syntheticWaveStack);
	    groundtruth.setStack(syntheticGroundTruthStack);
	    result.setSlice(t - 1);
	    groundtruth.setSlice(t - 1);
	    result.show();
	    groundtruth.show();	
	    result.updateAndDraw();
	    groundtruth.updateAndDraw();
	}

    }

    private ImageProcessor maskToProcessor(boolean[][] mask) {
	int width = mask.length;
	int height = mask[0].length;
	float[] newpixels = new float[width * height];
	FloatProcessor fp = new FloatProcessor(width, height);
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {

		newpixels[y * width + x] = mask[x][y] ? (float) normalInner
			.nextDouble() : (float) normalOuter.nextDouble();
	    }
	}
	fp.setPixels(newpixels);
	return fp;
    }

    private ImageProcessor maskToByteProcessor(boolean[][] mask) {
	int width = mask.length;
	int height = mask[0].length;
	byte[] newpixels = new byte[width * height];
	ByteProcessor fp = new ByteProcessor(width, height);
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {

		newpixels[y * width + x] = (byte) (mask[x][y] ? 0 : 255);
	    }
	}
	fp.setPixels(newpixels);
	return fp;
    }

    public int setup(String arg, ImagePlus imp) {
	this.imp = imp;
	return DOES_32 + DOES_8G + DOES_16 + NO_CHANGES + ROI_REQUIRED;

    }

}
