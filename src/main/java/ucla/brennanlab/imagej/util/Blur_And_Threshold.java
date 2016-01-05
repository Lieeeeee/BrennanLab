package ucla.brennanlab.imagej.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.StackStatistics;

public class Blur_And_Threshold implements PlugInFilter {
	ImagePlus imp;
	ImagePlus resultBlur;
	ImagePlus thresholdedBlur;
	int baseStart = 2, baseStop = 9;
	String arg;
	byte bitDepth;
	double sigma = -9;
	double blur = 9;

	private void Blur() {
		int nSlices = imp.getImageStackSize();

		ImageStack resultBlurStack = imp.createEmptyStack();
		ImageProcessor ip;
		GaussianBlur blurrer = new GaussianBlur();
		for (int slice = 1; slice <= nSlices; slice++) {
			ip = imp.getStack()
					.getProcessor(slice <= nSlices ? slice : nSlices)
					.duplicate();
			blurrer.blur(ip, blur);
			resultBlurStack.addSlice("", ip);
		}
		resultBlur = new ImagePlus("", resultBlurStack);
		ImageConverter conv32 = new ImageConverter(resultBlur);
		conv32.convertToGray32();

		return;
	}

	protected GenericDialog OptionSelector() {
		GenericDialog gd = new GenericDialog("Options", IJ.getInstance());
		gd.addMessage("Select from the options below:");
		gd.addNumericField("Baseline start slice:", baseStart, 0);
		gd.addNumericField("Stop slice:", baseStop, 0);
		gd.addNumericField("Threshold sigma:", sigma, 0);
		gd.addNumericField("Gaussian Blur sigma:", blur, 0);

		return (gd);

	}

	public void run(ImageProcessor ip) {
		GenericDialog gd = OptionSelector();
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		baseStart = (int) gd.getNextNumber();
		baseStop = (int) gd.getNextNumber();
		sigma = gd.getNextNumber();
		blur = gd.getNextNumber();
		Blur();
		Threshold();
		thresholdedBlur.show();

		return;

	}

	public int setup(String arg, ImagePlus imp) {
		this.arg = arg;
		this.imp = imp;
		baseStop = Math.min(imp.getImageStackSize(), 9);
		return STACK_REQUIRED + DOES_32;
	}

	private void Threshold() {
		ByteProcessor byteip = new ByteProcessor(resultBlur.getWidth(),
				resultBlur.getHeight());
		ImageStack bytestack = new ImageStack(resultBlur.getWidth(), resultBlur
				.getHeight());
		int nSlices = resultBlur.getImageStackSize();
		int dimension = resultBlur.getWidth() * resultBlur.getHeight();
		float[] pixels;
		byte[] newpixels;
		ImageStack baseBlurStack = resultBlur.createEmptyStack();
		for (int i = baseStart; i <= baseStop; i++) {

			ImageProcessor ip1 = resultBlur.getStack().getProcessor(
					i <= baseStop ? i : baseStop);
			ip1 = ip1.crop();
			baseBlurStack.addSlice("", ip1);

		}
		ImagePlus baseBlur = new ImagePlus("base", baseBlurStack);
		StackStatistics stackstats = new StackStatistics(baseBlur);
		double thresholdval = (stackstats.mean + sigma
				* stackstats.stdDev);
		for (int slice = 1; slice <= nSlices; slice++) {
			pixels = new float[dimension];
			newpixels = new byte[dimension];
			pixels = (float[]) resultBlur.getStack().getPixels(
					slice <= nSlices ? slice : nSlices);
			for (int j = 0; j < dimension; j++) {
				if (sigma < 0) {
					if ((pixels[j]) < thresholdval)
						newpixels[j] = (byte) (0 & 0xff);
					else
						newpixels[j] = (byte) (255 & 0xff);
				} else {
					if ((pixels[j]) > thresholdval)
						newpixels[j] = (byte) (0 & 0xff);
					else
						newpixels[j] = (byte) (255 & 0xff);
				}
			}
			byteip.setPixels(newpixels);
			bytestack.addSlice(slice + "", byteip);
		}

		thresholdedBlur = new ImagePlus("Thresholded", bytestack);
		return;
	}

}
