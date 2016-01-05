package ucla.brennanlab.imagej.util;

// Use stack statistics to find a good threshold value

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.StackStatistics;

public class Threshold_vs_Baseline implements PlugInFilter {
	ImagePlus imp;
	String arg;
	int baseStart = 2, baseStop = 9;
	int targetstart = 12, targetstop = 14;
	double blur = 9;
	double sigma, soffset = .3;

	protected GenericDialog OptionSelector() {
		GenericDialog gd = new GenericDialog("Options", IJ.getInstance());
		gd
				.addMessage("This plugin picks a threshold based upon the standard deviation of the baseline,\n and the number of standard deviations out the minimum of the target is");
		gd.addMessage("Select from the options below:");
		gd.addNumericField("Baseline start slice:", baseStart, 0);
		gd.addNumericField("Stop slice:", baseStop, 0);
		gd.addNumericField(
				"proportion of max standardized min-mean to threshold at",
				soffset, 2);
		gd.addNumericField("Gaussian Blur sigma:", blur, 1);
		gd.addNumericField("Target_start", targetstart, 0);
		gd.addNumericField("Target_stop", targetstop, 0);
		return (gd);

	}

	public void run(ImageProcessor ip) {
		GenericDialog gd = OptionSelector();
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		baseStart = (int) gd.getNextNumber();
		baseStop = (int) gd.getNextNumber();
		soffset = gd.getNextNumber();
		blur = gd.getNextNumber();
		targetstart = (int) gd.getNextNumber();
		targetstop = (int) gd.getNextNumber();

		// do our blurring first
		GaussianBlur blurrer = new GaussianBlur();
		ImageStack blurBaseStack = imp.createEmptyStack();
		ImageStack blurTargetStack = imp.createEmptyStack();
		ImageProcessor ip1;
		int nSlices = imp.getImageStackSize();

		ImageStack blurStack = imp.createEmptyStack();

		for (int slice = 1; slice <= nSlices; slice++) {
			ip1 = imp.getStack().getProcessor(
					slice <= nSlices ? slice : nSlices).duplicate();
			blurrer.blur(ip1, blur);
			blurStack.addSlice("Slice " + slice, ip1);
			if (slice >= baseStart && slice <= baseStop)
				blurBaseStack.addSlice("Slice " + slice, ip1);
			if (slice >= targetstart && slice <= targetstop)
				blurTargetStack.addSlice("Slice " + slice, ip1);
			IJ.showProgress(slice, nSlices);

		}

		ImagePlus blurBase = new ImagePlus("", blurBaseStack);
		ImagePlus blurTarget = new ImagePlus("", blurTargetStack);

		StackStatistics basestats = new StackStatistics(blurBase);
		StackStatistics targetstats = new StackStatistics(blurTarget);

		double sdrange = (targetstats.min - basestats.mean) / basestats.stdDev;

		IJ.showStatus("standard min range: " + sdrange);

		double sigma = sdrange - sdrange * soffset;
		double threshold = basestats.mean + sigma * basestats.stdDev;

		ByteProcessor byteip = new ByteProcessor(imp.getWidth(), imp
				.getHeight());
		ImageStack bytestack = new ImageStack(imp.getWidth(), imp.getHeight());

		float[] pixels;
		byte[] newpixels;
		int dimension = imp.getWidth() * imp.getHeight();

		for (int slice = 1; slice <= nSlices; slice++) {
			pixels = new float[dimension];
			newpixels = new byte[dimension];
			pixels = (float[]) blurStack.getPixels(slice <= nSlices ? slice
					: nSlices);
			for (int j = 0; j < dimension; j++) {
				if (sigma < 0) {
					if ((pixels[j]) < threshold)
						newpixels[j] = (byte) (0 & 0xff);
					else
						newpixels[j] = (byte) (255 & 0xff);
				} else {
					if ((pixels[j]) > threshold)
						newpixels[j] = (byte) (0 & 0xff);
					else
						newpixels[j] = (byte) (255 & 0xff);
				}
			}
			byteip.setPixels(newpixels);
			bytestack.addSlice(slice + "", byteip);
		}

		ImagePlus thresholdedBlur = new ImagePlus("Thresholded", bytestack);
		thresholdedBlur.show();

		return;

	}

	public int setup(String arg, ImagePlus imp) {
		this.arg = arg;
		this.imp = imp;
		baseStop = Math.min(imp.getImageStackSize(), 9);
		return STACK_REQUIRED + DOES_32;
	}

}