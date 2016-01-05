package ucla.brennanlab.imagej.csd;

/*
 The use of this plug-in requires user intervention. This plug-in will ask 
 the user to pick out a CSD start frame as well as an approximate location
 for the CSD start locus. The user will then be asked to pick a frame to
 stop analysis. This plugin will then interface with FeatureJ Derivatives
 to compute the time derivative of the image within the user defined frames 
 + some padding on either ends. The user's input constitutes the prior 
 information in this model.
 */
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.filter.RankFilters;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

public class CSD_peak_isolator implements PlugInFilter {
	public static double max(double[] t) {
		double maximum = t[0]; // start with the first value
		for (int i = 1; i < t.length; i++) {
			if (t[i] > maximum) {
				maximum = t[i]; // new maximum
			}
		}
		return maximum;
	}// end method max

	public static int max(short[] t) {
		int maximum = t[0]; // start with the first value
		for (int i = 1; i < t.length; i++) {
			if (t[i] > maximum) {
				maximum = t[i]; // new maximum
			}
		}
		return maximum;
	}// end method max

	public static double mean(byte[] p) {

		double sum = 0; // sum of all the elements
		for (int i = 0; i < p.length; i++) {
			sum += p[i] & 0xFF;
		}
		return sum / p.length;
	}// end method mean

	public static double min(double[] t) {
		double min = t[0]; // start with the first value
		for (int i = 1; i < t.length; i++) {
			if (t[i] < min) {
				min = t[i]; // new maximum
			}
		}
		return min;
	}// end method max

	public static int min(short[] t) {
		int min = t[0]; // start with the first value
		for (int i = 1; i < t.length; i++) {
			if (t[i] < min) {
				min = t[i]; // new maximum
			}
		}
		return min;
	}// end method max

	public static double sd(double[] vals) {
		double sum = 0;
		double sum2 = 0;
		for (int j = 0; j < vals.length; j++) {
			sum += vals[j];
			sum2 += Math.pow(vals[j], 2.0);
		}
		double mean = sum / vals.length;
		double var = sum2 / vals.length - Math.pow(mean, 2);
		return Math.pow(var, 0.5);
	}

	public static int sum(int[] t) {
		int sum = 0;
		for (int j = 0; j < t.length; j++) {
			sum += t[j];
		}
		return sum;
	}

	private ImagePlus imp; // input image
	public ImagePlus binaryWave;
	String arg;
	private int N;
	private int width;
	private int height;
	private int timePt; // user identified time point
	private int frameSkip = 1;
	public int arrivalTime;
	public double area;
	short[] leading; // slice of beginning of CSD
	short[] peak; // slice of peak
	short[] lagging;
	short[] sliceRange;

	short[] meanWavePeak;

	int keep = 100;

	double peakMagnitude;

	private int critical = 1024;

	private boolean regularize = true;

	private boolean showArrivalTimes = false;

	private boolean showMovie = true;

	public CSD_peak_isolator() {

	}

	public CSD_peak_isolator(ImagePlus imp, int skip, int keep) {
		this.imp = imp;
		this.N = imp.getStackSize();
		this.width = imp.getWidth();
		this.height = imp.getHeight();
		this.leading = new short[width * height];
		this.peak = new short[width * height];
		this.lagging = new short[width * height];
		this.sliceRange = new short[] { 0, 0, 0 };
		this.meanWavePeak = new short[] { 0, 0, 0 };
		this.frameSkip = skip;
		this.keep = keep;

	}

	public short[] findLargestPeak(double[] oisVals) {
		/*
		 * return a vector with 1) start of CSD brightening, 2) peak of
		 * brightening
		 */
		short peak;
		// First smooth the original values
		double[] smoothed = smooth(oisVals);

		// this is a lagging delta (how much the smoothed value at t is
		// different from that at t-5
		double[] deltas;
		deltas = new double[smoothed.length];
		for (int i = 0; i < 8; i++) { // don't let CSD occur in the first 10
			// frames
			deltas[i] = 0;
		}
		for (int i = 8; i < smoothed.length; i++) {
			deltas[i] = smoothed[i] - smoothed[i - 8];
		}

		// now find the index with the largest change within 10 frames.
		// use the user input if available
		short biggestchange;
		if (timePt < 1 || timePt > oisVals.length) {
			biggestchange = whichIsBiggest(deltas);
		} else
			biggestchange = (short) timePt;

		// now search in the vicinity of this index and find the real peak
		while (biggestchange <= smoothed.length - 2
				&& smoothed[biggestchange + 1] > smoothed[biggestchange]) {
			biggestchange++;
		}

		peak = biggestchange;
		// find start
		short start = peak;
		while (start > 0 && smoothed[start - 1] < smoothed[start]) {
			start--;
		}
		short end = peak;
		while (end < smoothed.length - 2 && smoothed[end + 1] < smoothed[end]) {
			end++;
		}
		// If we don't see CSD-like changes, return the last time point in the
		// stack.
		// this applies also if the pixel doesn't change much..
		// the change associated with CSD is big, so we expect it to be at least
		// as big as the range

		double sigma = sd(oisVals);
		if (this.sliceRange[1] != 0) {
			if ((sigma == 0)
					|| (oisVals[peak] - oisVals[start]) < 0.20 * this.peakMagnitude
					|| Math.abs(peak - start) < 10) {
				return new short[] { (short) oisVals.length,
						(short) oisVals.length, (short) oisVals.length };
			}
		}
		return new short[] { (short) (start + 1), (short) (peak + 1),
				(short) (end + 1) };
	}

	public ImagePlus getLeadingTimes() {

		//
		ImageProcessor leadingLevelSetIP = new ShortProcessor(width, height,
				this.leading, null);
		if (this.regularize) {
			RankFilters filter = new RankFilters();
			filter.rank(leadingLevelSetIP, 1, RankFilters.MEDIAN);
		}
		//

		ImagePlus leading = new ImagePlus("CrossingTimes", leadingLevelSetIP);
		return leading;
	}

	public void isolatePeaks() {
		/* First take z projection mean for entire stack, then find csd */
		ImageStack stack = imp.getStack();
		double[] meanByFrame;
		meanByFrame = new double[N];
		byte[] pixels;
		for (int i = 0; i < N; i++) {
			pixels = (byte[]) stack.getPixels(i + 1);
			meanByFrame[i] = mean(pixels);
		}
		this.meanWavePeak = findLargestPeak(meanByFrame);

		this.peakMagnitude = (meanByFrame[meanWavePeak[1] - 1] - meanByFrame[meanWavePeak[0] - 1]);
		this.sliceRange[0] = (short) Math.max(1, this.meanWavePeak[0] - 50);
		this.sliceRange[1] = this.meanWavePeak[1];
		this.sliceRange[2] = (short) Math.min(N, this.meanWavePeak[2] + 30);

		double[] fixedPtVals;
		short[] locPeak;

		for (int h = 0; h < height; h++) {
			IJ.showStatus("Detecting pointwise peaks...");
			IJ.showProgress(h + 1, height);
			for (int w = 0; w < width; w++) {
				fixedPtVals = new double[this.sliceRange[2]
						- this.sliceRange[0] + 1]; // the pixel values at this
				// fixed point
				locPeak = new short[3];

				for (int slice = this.sliceRange[0]; slice <= this.sliceRange[2]; slice++) {
					fixedPtVals[slice - this.sliceRange[0]] = ((byte[]) stack
							.getPixels(slice))[h * width + w] & 0xFF;
				}
				locPeak = findLargestPeak(fixedPtVals);
				this.leading[h * width + w] = (short) (locPeak[0] + this.sliceRange[0]);
				this.peak[h * width + w] = (short) (locPeak[1] + this.sliceRange[0]);
				this.lagging[h * width + w] = (short) (locPeak[2] + this.sliceRange[0]);
			}
		}
	}

	public ImagePlus makeBinaryWave() {
		// This takes this.leading and makes it into a progressive threshold
		// binary stack

		ImageStack LeadingBinary = imp.createEmptyStack();
		ByteProcessor currentSliceIP;
		byte[] pixels;
		ImageProcessor leadingLevelSetIP = new ShortProcessor(width, height,
				this.leading, null);
		ImageProcessor laggingLevelSetIP = new ShortProcessor(width, height,
				this.lagging, null);
		ImageProcessor peakLevelSetIP = new ShortProcessor(width, height,
				this.peak, null);
		RankFilters filter = new RankFilters();

		if (regularize) {
			filter.rank(leadingLevelSetIP, 1, RankFilters.MEDIAN);
			filter.rank(laggingLevelSetIP, 1, RankFilters.MEDIAN);
			filter.rank(peakLevelSetIP, 4, RankFilters.MEDIAN);
		}
		short[] laggingpx = (short[]) laggingLevelSetIP.getPixels();
		short[] leadingpx2 = (short[]) leadingLevelSetIP.getPixels();
		short[] peakpx = (short[]) peakLevelSetIP.getPixels();

		boolean commenced = false;

		short[] pixels2 = new short[leadingpx2.length];
		int firsttime = min(leadingpx2);
		int lasttime = max(this.leading);

		while (firsttime < lasttime && !commenced) { // find the first slice of
			// interest
			int count = 0;
			IJ.showStatus("Checking slice " + firsttime + " for CSD...");

			for (int l = 0; l < leadingpx2.length; l++) {
				if (leadingpx2[l] <= firsttime) {
					pixels2[l] = 0;
					count++;
				} else
					pixels2[l] = 255;
			}

			/*
			 * leadingSetIP=new ShortProcessor(width,height,this.leading,null);
			 * filter.rank(leadingSetIP, 1, RankFilters.MEDIAN);
			 * 
			 * short[] waveindicator = (short[]) leadingSetIP.getPixels();
			 * 
			 * for(int m=0;m<waveindicator.length;m++){ //
			 * if(waveindicator[m]==0) count++; }
			 */

			if (count < this.critical)
				firsttime++;
			else {
				commenced = true;
				IJ.showStatus("CSD found");
			}

		}

		int kept = 0;
		this.arrivalTime = firsttime;
		int areacnt = 0;

		for (int i = firsttime; i < lasttime; i++) {
			if (kept >= keep)
				break;
			if ((i - firsttime) % this.frameSkip == 0) {

				if (regularize)
					IJ.showStatus("Regularizing...");
				else
					IJ.showStatus("Assembling stack...");
				IJ.showProgress(i - firsttime, firsttime + lasttime);
				pixels = new byte[leadingpx2.length];
				for (int j = 0; j < laggingpx.length; j++) {
					if (i >= laggingpx[j]) {
						pixels[j] = (byte) 220; // past acute

					} else if (i >= peakpx[j])
						pixels[j] = (byte) 80; // acute recovery
					else if (leadingpx2[j] <= i) {
						pixels[j] = (byte) 0; // acute brightening 220, 80, 0
					} else {
						pixels[j] = (byte) 255; // pre-acute
						if (i > lasttime - frameSkip - 1 || kept == keep - 1)
							areacnt++;
					}
				}

				// rank filter

				currentSliceIP = new ByteProcessor(width, height, pixels, null);
				if (this.regularize)
					filter.rank(currentSliceIP, 1, RankFilters.MEDIAN);
				LeadingBinary.addSlice(imp.getStack().getSliceLabel(i)
						+ " time frame " + i, currentSliceIP);
				kept++;
			}
			this.area = 1 - (double) areacnt / width / height;

		}
		ImagePlus wave = new ImagePlus("CSD Wave", LeadingBinary);
		return (wave);

	}

	public void run(ImageProcessor ip) {

		isolatePeaks();
		this.binaryWave = makeBinaryWave();
		if (showMovie) {
			binaryWave.show();
			binaryWave.updateAndDraw();
			binaryWave.unlock();
		}
		if (showArrivalTimes) {
			ImagePlus arrivalTimes = getLeadingTimes();
			arrivalTimes.show();
			arrivalTimes.updateAndDraw();
			arrivalTimes.unlock();
		}
		this.imp.unlock();

	}

	public void runWithTxy(ImageProcessor ip) {
		this.showArrivalTimes = true;
		run(ip);
	}

	public void setRegularization(boolean regularize) {
		this.regularize = regularize;
	}

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		this.arg = arg;
		this.N = imp.getStackSize();
		this.width = imp.getWidth();
		this.height = imp.getHeight();
		this.leading = new short[width * height];
		this.peak = new short[width * height];
		this.lagging = new short[width * height];
		this.sliceRange = new short[] { 0, 0, 0 };
		this.meanWavePeak = new short[] { 0, 0, 0 };

		if (arg.length() > 0) {
			String[] args = arg.split(";");
			for (int i = 0; i < args.length; i++) {
				String[] param = args[i].split("=");
				if (param[0].equals("timePt"))
					timePt = new Integer(param[1]);
				else if (param[0].equals("frameSkip"))
					frameSkip = new Integer(param[1]);
				else if (param[0].equals("showArrivalTimes"))
					showArrivalTimes = new Boolean(param[1]);
				else if (param[0].equals("showMovie"))
					showMovie = new Boolean(param[1]);
				else if (param[0].equals("regularize"))
					regularize = new Boolean(param[1]);

			}
		}

		return STACK_REQUIRED + NO_CHANGES + DOES_8G;
	}

	public double[] smooth(double[] data) {
		// savitzky-golay coefficients
		int[] sgcoefficients = new int[] { -21, 14, 39, 54, 59, 54, 39, 14, -21 };
		double sgnorm = sum(sgcoefficients);
		int width = (int) Math.floor(sgcoefficients.length / 2);
		double smoothed[] = new double[data.length];
		for (int j = 0; j < data.length; j++) {
			smoothed[j] = 0;
			for (int k = 0; k < sgcoefficients.length; k++) {
				smoothed[j] += sgcoefficients[k]
						* data[Math.min(data.length - 1, Math.max(0, j
								- width + k))] / sgnorm;
			}
		}
		return smoothed;
	}

	public short whichIsBiggest(double[] vals) {
		int largestndx = 0;
		for (int j = 1; j < vals.length; j++) {
			if (vals[j] > vals[largestndx]) {
				largestndx = j;
			}
		}
		return (short) largestndx;
	}
}
