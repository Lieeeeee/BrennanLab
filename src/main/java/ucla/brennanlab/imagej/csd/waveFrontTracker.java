package ucla.brennanlab.imagej.csd;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.Roi;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.Color;
import java.awt.Frame;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;
import ucla.brennanlab.imagej.graphcuts.GraphCutSegmenter;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;
import ucla.brennanlab.imagej.util.levelsets.ShapePrior;
import ucla.brennanlab.imagej.util.stochastic.IntensityModel;
import ucla.brennanlab.imagej.util.stochastic.LaplaceMixture;
import ucla.brennanlab.imagej.util.stochastic.SpeedField;

/**
 * Boundary segmentation using the kriging speed model
 * 
 * @author Josh Chang
 * 
 */
public class waveFrontTracker implements PlugInFilter {
    ImagePlus imp;
    int timesteps;
    int slices;
    int firstslice;
    int gridwidth = 5;// reduce the dimensions down to 2x2 blocks
    int gridheight = 5;
    int bandwidth = 30; // band for interpoldation, all outside is set to

    double initialLengthPenalty = 20;
    double lengthPenalty = 12;
    double alpha = 2;
    boolean showProgress = true;
    boolean showPredictions = false;
    boolean colorResult = true;
    boolean useManualInitializationPrior = false;
    Roi userDrawnRoi;
    /**
     * the number of importance samples is J*(2L+1). L is usually the dimension
     * of the stochastic state, here it will be the number of points at which we
     * will estimate the speed, interpolating speeds using these points
     */
    static final int numContourSamplesJ = 2;
    int speedSamples = 16;
    double conversion = 10, freq = 1, priorSpeedMean = 3, priorSpeedSD = 1;
    boolean verbose = true;
    int width;
    int height;
    int innerBandWidth = 25;
    SpeedField maxAposteriorSpd; // filtered speed field

    final double lambda = 0.1;
    ArrayList<Roi> roisOfSegmentations;
    boolean csdHasStarted = false;
    Color roiColor = new Color(255, 105, 180,255); // Color to draw the ROI

    public void run(ImageProcessor ip) {
	
	IJ.log("Starting CSD Wavefront segmentation on " + imp.getTitle());
	Frame logFrame = WindowManager.getFrame("Log");
	if (logFrame != null) {
	    logFrame.setSize(1024, 480);
	}

	// Obtain user input
	GenericDialog gd = pluginDialog();
	gd.showDialog();
	if (gd.wasCanceled()) {
	    return;
	}

	conversion = gd.getNextNumber();
	freq = gd.getNextNumber();
	double inputspeed = gd.getNextNumber();
	double inputstder = gd.getNextNumber();
	alpha = gd.getNextNumber();
	lengthPenalty = gd.getNextNumber();
	initialLengthPenalty = gd.getNextNumber();
	priorSpeedMean = 100 * inputspeed / 6 / conversion / freq;
	priorSpeedSD = 100 * inputstder / 6 / conversion / freq;
	this.speedSamples = (int) gd.getNextNumber();
	this.innerBandWidth = (int) gd.getNextNumber();
	this.showPredictions = gd.getNextBoolean();
	
	if(userDrawnRoi !=null)
	    this.useManualInitializationPrior = gd.getNextBoolean();

	ImageWindow iw = WindowManager.getCurrentWindow();
	
	/**********************************
	 * keep user clicks // on other images // from screwing // things up
	 **********************************/
	int currentSlice = imp.getCurrentSlice(); // the current slice
	roisOfSegmentations = new ArrayList<Roi>(imp.getSlice() + 100 < imp
		.getStackSize() ? imp.getSlice() + 100 : imp.getStackSize());

	IJ.log("Running segmentation with params: nu_0 "
		+ this.initialLengthPenalty + " nu: " + this.lengthPenalty
		+ " alpha: " + this.alpha);
	imp.unlock();

	/*********************************************************************
	 * Initialize and display the ROI manager and log viewer
	 *********************************************************************/
	RoiManager rm = RoiManager.getInstance();
	if (rm != null)
	    rm.setVisible(true);
	else
	    rm = new RoiManager();

	GraphCutSegmenter gcs = new GraphCutSegmenter(ip);

	ImageStack gSegmentStack = imp.createEmptyStack();
	ImagePlus gSegment = imp.createImagePlus();
	gSegment.setTitle(imp.getTitle() + "_segmentation");

	this.maxAposteriorSpd.setPriorMeanSpeed((float) this.priorSpeedMean);
	this.maxAposteriorSpd.setPriorVarSpeed((float) this.priorSpeedSD);

	/*****************************************************************
	 * Start processing totalSlices here.
	 * 
	 * The following is code for a null shape:
	 * 
	 * ShapePrior nullShapeKernelPrior = new
	 * ShapePrior( width, height);
	 * 
	 *ImplicitShape2D nullShape = nullShapeKernelPrior.shapes.get(0);
	 * 
	 *****************************************************************/

	int firstCSDslice = currentSlice;
	int lastCSDslice = currentSlice + 100 < imp.getStackSize() ? currentSlice + 100
		: imp.getStackSize();
	RandomEngine engine = new cern.jet.random.engine.MersenneTwister(
		new java.util.Date());
	final Normal normalRV = new Normal(this.priorSpeedMean,
		this.priorSpeedSD, engine);

	NumberFormat fmt = new DecimalFormat("##");

	while (currentSlice <= lastCSDslice) {
	    if (IJ.escapePressed()) {
		IJ.beep();
		IJ.log("Aborted");
		return;
	    }
	    WindowManager.setCurrentWindow(iw); // In case the user clicked on
	    // another window
	    imp.setSlice(Math.max(1, currentSlice));
	    ImageProcessor currentImageProcessor = imp.getStack().getProcessor(
		    currentSlice <= lastCSDslice ? currentSlice : lastCSDslice);
	    currentImageProcessor = currentImageProcessor.crop();

	    /*************************************************************************************
	     * If we determine that CSD has started, run this conditional
	     *************************************************************************************/

	    if (this.useManualInitializationPrior) {
		// Adjust the length penalty perhaps as a function of time?
		this.lengthPenalty = 4;
	    }

	    if (csdHasStarted) {

		ImplicitShape2D prevLS =new ImplicitShape2D(gcs.returnMask());
		if (gcs.fractionInner() < 10.0 / width / height) {
		    //prevLS = 
		}

		/****************************************************************
		 * Time to propagate through the level-set equation using fast
		 * marching. from this calculation, determine the signed
		 * distance functions, then reconstruct the prior using kernel
		 * density estimation, does multithreaded soooon :)
		 ***************************************************************/

		ArrayList<ImplicitShape2D> positions = new ArrayList<ImplicitShape2D>(
			speedSamples);
		Roi[] prevPos = new Roi[roisOfSegmentations.size()];
		for (int p = 0; p < prevPos.length; p++) {
		    prevPos[p] = (Roi) roisOfSegmentations.get(p).clone();

		}

		ImplicitShape2D meanNext = new ImplicitShape2D(prevLS
				.getMask());
		meanNext.advectUniformSpeed(this.priorSpeedMean, 1 / this.freq);
		
		// advect nonuniform speed damn it!!

		Roi nextRoi = meanNext
			.getRoi(false);
		imp.setRoi(nextRoi);

		nextRoi.setName("Predicted mean position " + currentSlice);
		rm.addRoi(nextRoi);
		
		SpeedField sfield = new SpeedField(this.width,this.height);

		double[] effectiveWeights = new double[this.speedSamples];
		double speed;
		for (int i = 0; i < this.speedSamples; i++) {
		    speed = normalRV.nextDouble(priorSpeedMean, priorSpeedSD);

		    positions.add(new ImplicitShape2D(
			   prevLS.getMask()));
		    effectiveWeights[i] = 1;// redo the weighting soon
		    
		    positions.get(i).advectUniformSpeed(speed, 1 / this.freq);

		    Roi pred = (Roi) positions.get(i).getRoi(true);

		    if (pred == null)
			continue;

		    if (showPredictions) {
			pred
				.setName("Sample " + (i + 1)
					+ " of predicted position of "
					+ (currentSlice));
			rm.addRoi(pred);
		    }

		}

		effectiveWeights = normalize(effectiveWeights);

		/*************************************************************
		 * Kernel density estimation on the shapes to determine prior
		 *************************************************************/

		ShapePrior skde = new ShapePrior(
			positions, effectiveWeights);

		//skde.setAlpha(alpha);
		IJ.log("Prediction distance variance " + skde.tauSquared);

		IJ.log("... done computing priors");
		IJ.log("... computing maxflow for slice " + currentSlice
			+ "...");
		long startTime = System.nanoTime();
		float minE;
		float tempE;
		long endTime;
		IntensityModel s = new LaplaceMixture();
		s.Infer(currentImageProcessor, meanNext
			.getMask());

		gcs = new GraphCutSegmenter(currentImageProcessor);
		gcs.setLengthPenalty((float)this.lengthPenalty);
		
		/********************************************************************************
		 * Segment NOW
		 *********************************************************************************/
		
		
		try {
		    gcs.setNodeWeights(skde, meanNext, s);
		    gcs.setEdgeWeights(skde, meanNext);
		    minE = gcs.relaxEnergy();
		    /**
		     * MM iterations
		     */
		    s.Infer(currentImageProcessor, gcs
			    .returnMask());
		    gcs.setEdgeWeights(skde, gcs.getLevelSet());
		    long currentTime = startTime;
		    int i = 1;
		    /*
		     * Every 5 iterations or so, let's show the progress!!
		     */

		    while (true) {
			if (IJ.escapePressed()) {
			    IJ.beep();
			    IJ.log("Aborted");
			    return;
			}
			gcs.gc.reset();
			gcs.setNodeWeights(skde, gcs.getLevelSet(), s);
			gcs.setEdgeWeights(skde, gcs.getLevelSet());
			tempE = gcs.relaxEnergy();
			currentTime = System.nanoTime();
			IJ.log("\t     Iter " + i + " Elapsed time: "
				+ ((float) (currentTime - startTime) / 1000000)
				+ " ms" + " Energy " + tempE);
			if (minE - tempE < 1 && i > 2 || Double.isNaN(minE)
				|| Double.isInfinite(minE)) {
			    minE = tempE;
			    if (Double.isNaN(minE) || Double.isInfinite(minE)) {
				IJ.log("Encountered invalid energy");
			    }

			    break;
			}
			minE = tempE;
			// if (i % 3 == 0) {
			try {
			    Roi currentProgress = (Roi) gcs.getLevelSet().getRoi(
				    false).clone();
			    currentProgress.setStrokeColor(roiColor);

			    currentProgress.setName("Position_" + currentSlice
				    + "_iteration_" + i);
			    imp.setRoi(currentProgress);
			    rm.addRoi(currentProgress);

			} catch (NullPointerException e) {
			    csdHasStarted = false;

			}

			// }
			i++;
		    }

		} finally {
		    endTime = System.nanoTime();
		}
		Roi segmentedRoi = gcs.returnRoi();
		segmentedRoi.setStrokeColor(roiColor);

		if (segmentedRoi == null) {
		    segmentedRoi = meanNext
			    .getRoi(false);
		}
		if (segmentedRoi != null) {
		    if (currentSlice == firstCSDslice + 1) { // Draw first slice
			gSegmentStack.addSlice("Slice_" + firstCSDslice, imp
				.getStack().getProcessor(firstCSDslice)
				.duplicate().convertToRGB());
			roisOfSegmentations.get(0).setStrokeColor(this.roiColor);
			roisOfSegmentations.get(0).setStrokeWidth(3);
			gSegmentStack.getProcessor(1).setColor(this.roiColor);
			gSegmentStack.getProcessor(1).drawRoi(
				roisOfSegmentations.get(0));
		    }
		    // Draw other slices
		    gSegmentStack.addSlice("Slice_" + currentSlice,
			    currentImageProcessor.duplicate().convertToRGB());

		    gSegment.setStack(gSegmentStack);
		    gSegment.show();
		    gSegment.setSlice(gSegment.getStackSize());
		    gSegment.updateAndDraw();

		    roisOfSegmentations.add((Roi) segmentedRoi.clone());
		    roisOfSegmentations.get(currentSlice - firstCSDslice)
			    .setName("Position_" + currentSlice);
		    roisOfSegmentations.get(currentSlice - firstCSDslice)
			    .setStrokeColor(this.roiColor);
		    roisOfSegmentations.get(currentSlice - firstCSDslice)
			    .setStrokeWidth(3);
		
		    rm.addRoi(roisOfSegmentations.get(currentSlice
			    - firstCSDslice));


		    gSegmentStack.getProcessor(currentSlice-firstCSDslice+1).setColor(this.roiColor);
		   ((ColorProcessor)  gSegmentStack
			    .getProcessor(currentSlice - firstCSDslice + 1))
			    .drawRoi(
				    roisOfSegmentations.get(currentSlice
					    - firstCSDslice));
		    gSegment.setSlice(currentSlice - firstCSDslice + 1);

		} else {
		    IJ
			    .log("Something went wrong, I was unable to obtain a segmentation...");

		    lastCSDslice = currentSlice - 1;
		    continue;
		}

		long duration = endTime - startTime;
		double dist = skde.computeDistance(gcs.getLevelSet(), meanNext)
			*skde.tauSquared;
		IJ.log("Computed maxflow energy of " + minE + " in "
			+ (float) duration / 1000000 + " ms");
		IJ
			.log("Normalized distance between mean prediction and result: "
				+ dist);

		if (Double.isNaN(minE) || minE == 0 || Double.isNaN(dist)) {
		    lastCSDslice = currentSlice - 1;
		    continue;
		}

		ImplicitShape2D currentPosition = gcs.getLevelSet();

		// perhaps do some sanity checks
		
		this.maxAposteriorSpd.addArrival(currentPosition);

		this.maxAposteriorSpd.updateGMRF();
		
		this.priorSpeedMean = 0.9 * this.maxAposteriorSpd.latestSpeed
			+ 0.1 * this.priorSpeedMean;
		if (1/s.getPosteriorPrecision(false) > s.getPosteriorMean(true)) {
		    IJ
			    .log("\t\t Excessive variability in the image detected, may be movement");
		     this.priorSpeedMean = this.priorSpeedMean * 0.8;
		}
		this.priorSpeedSD = this.priorSpeedMean * 0.5; // / 1.25;

		IJ.log("Current average wavespeed is "
			+ this.maxAposteriorSpd.latestSpeed / 100 * 6
			* conversion * freq + " mm/s"); // current wavespeed

	    } else {
		/****************************************************
		 * CSD maybe hasn't started yet
		 *****************************************************/

		if (this.useManualInitializationPrior) {
		    this.userDrawnRoi.setName("User_drawn_template_slice_"
			    + currentSlice);
		    rm.addRoi(this.userDrawnRoi);
		    ImplicitShape2D manualLS = (new ImplicitShape2D(
			    this.userDrawnRoi, width, height));
		    ArrayList<ImplicitShape2D> userDefinedShapes = new ArrayList<ImplicitShape2D>(
			    1);
		    userDefinedShapes.add(manualLS);
		    double[] weights = new double[1];
		    weights[0] = 1;
		    ShapePrior userKDE = new ShapePrior(
			    userDefinedShapes, weights);
		    userKDE.setlambda(2.0);
		    userKDE.tauSquared = 1;
		    IntensityModel s = new LaplaceMixture();
		    s.Infer(currentImageProcessor, manualLS.getMask());

		    gcs = new GraphCutSegmenter(currentImageProcessor, 1);

		    gcs.setIntensityModel(s);
		    gcs.setEdgeWeights(userKDE, manualLS);
		    gcs.setNodeWeights(userKDE, manualLS, s);
		    double minE = gcs.relaxEnergy();
		    IJ.log(minE + " " + gcs.fractionInner());
		    s.Infer(currentImageProcessor, gcs.returnMask());

		} else {
		    IJ.log("... computing naive maxflow for slice "
			    + currentSlice + " ... ");
		    long startTime = System.nanoTime();
		    long endTime;
		    try {
			gcs = this.naivelyFindSegmentation(
				currentImageProcessor, 3);

		    } finally {
			endTime = System.nanoTime();
		    }

		    long duration = endTime - startTime;

		    IJ.log("Computed maxflow energy of " + gcs.Energy + " in "
			    + (float) duration / 1000000 + " ms");
		}
		// double inner = gcs.getS().innerarea;
		if (gcs.fractionInner() > 0.01) {
		    /**
		     * Provisionally say that the CSD has been initialized, but
		     * we might get kicked back here again if the subsequent
		     * frames yield inadmissible wave speed profiles
		     */
		    csdHasStarted = true;
		    IJ.log("Found a " + fmt.format(gcs.fractionInner() * 100)
			    + "% blob that might be the start of CSD at slice "
			    + currentSlice + "...");
		    Roi initialROI = gcs.returnRoi();
		    initialROI.setName("Position_" + currentSlice);
		    rm.addRoi(initialROI);
		    firstCSDslice = currentSlice;
		    roisOfSegmentations.add((Roi) initialROI.clone());
		    this.maxAposteriorSpd.addArrival(gcs.getLevelSet());

		} else {
		    csdHasStarted = false;
		    IJ.log("Found no trace of CSD, skipping to next slice\n");
		}

	    }

	    currentSlice++;

	}

	IJ.setSlice(currentSlice);

	gSegment.setOverlay(null);
	imp.setOverlay(null);

	IJ.log("Done segmenting, check out the output!");

	IJ
		.log("Computing arrival times and speed map for the filtered trajectory...");
	float[][] arrivalTimes = this.maxAposteriorSpd.getArrivalTimes();
	//GaussianBlur gb = new GaussianBlur();

	ImageProcessor timeIP = new FloatProcessor(arrivalTimes);

	(new ImagePlus("Arrival_times", timeIP)).show();
	arrivalTimes = new float[width][height];
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		arrivalTimes[x][y] = timeIP.getPixelValue(x, y);
	    }
	}

	float[][] speeds = this.maxAposteriorSpd.arrivalsToSpeed(arrivalTimes);
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		speeds[x][y] = (float) (speeds[x][y] / 100 * 6 * conversion * freq);
	    }
	}
	ImageProcessor speedIP = new FloatProcessor(speeds);
	speedIP.setMinAndMax(this.priorSpeedMean-this.priorSpeedSD/2, this.priorSpeedMean * 1.5);
	ImagePlus speedImage = new ImagePlus("Speed_field_mm_per_min", speedIP);
	//gb.blur(speedIP, 1);
	speedImage.show();

    }

    private double[] normalize(double[] values) {
	double sum = 0;
	for (int i = 0; i < values.length; i++) {
	    sum += values[i];
	}
	for (int i = 0; i < values.length; i++) {
	    values[i] /= sum;
	}
	return values;
    }

    public GraphCutSegmenter naivelyFindSegmentation(ImageProcessor ip,
	    int iterations) {
	IntensityModel s = new LaplaceMixture();
	GraphCutSegmenter gcs = new GraphCutSegmenter(ip);
	gcs.setIntensityModel(s);
	boolean[][] mask = new boolean[ip.getWidth()][ip.getHeight()];
	ImageProcessor ip2 = ip.duplicate();
	GaussianBlur gb = new GaussianBlur();
	gb.blur(ip2, 2);
	ImageStatistics is = ImageStatistics.getStatistics(ip2,
		ij.measure.Measurements.MEAN + ij.measure.Measurements.STD_DEV,
		null);

	double mean = is.mean;
	double sd = is.stdDev;

	// Just some initial values
	for (int x = 0; x < ip.getWidth(); x++) {
	    for (int y = 0; y < ip.getHeight(); y++) {
		mask[x][y] = ip2.getPixelValue(x, y) > mean + sd ? true : false;
	    }
	}

	s.Infer(ip2, mask);
	gcs.setLengthPenalty((float) this.initialLengthPenalty);
	int i = 1;
	gcs.setEdgeWeights();
	while (i++ < iterations) {
	    gcs.setNodeWeights(s);
	    IJ.log("Energy " + gcs.relaxEnergy());
	    s.Infer(ip, gcs.returnMask());
	}
	return gcs;

    }

    public int setup(String arg, ImagePlus imp) {
	this.imp = imp;
	this.timesteps = imp.getStackSize();
	this.slices = imp.getStackSize();
	this.width = imp.getWidth();
	this.height = imp.getHeight();
	this.maxAposteriorSpd = new SpeedField(this.width, this.height);
	this.userDrawnRoi = imp.getRoi();

	return DOES_16 + DOES_8G + DOES_32 + NO_CHANGES;
    }

    protected GenericDialog pluginDialog() {
	GenericDialog gd = new GenericDialog("Segmentation Options", IJ
		.getInstance());
	gd.addNumericField("Microns per pixel", 10, 3);
	gd.addNumericField("Acquisition frequency", 1, 1, 3, "Hz");
	gd.addNumericField("Prior mean wave speed", this.priorSpeedMean, 2, 5,
		"mm per minute");
	gd.addNumericField("standard error", this.priorSpeedSD, 2, 5,
		"mm per minute");
	gd.addNumericField("Shape penalty exponent", this.alpha, 0);
	gd.addNumericField("Length penalty", this.lengthPenalty, 0);
	gd.addNumericField("Initial Length penalty", this.initialLengthPenalty,
		0);
	gd.addNumericField("Number of predictions", this.speedSamples, 0);
	gd.addNumericField("Width of inner region for stats",
		this.innerBandWidth, 0);
	gd.addCheckbox("Show predicted positions", false);

	if (userDrawnRoi != null) {
	    gd.addCheckbox("Use userdrawn ROI", true);
	    gd
		    .addMessage("Note: User-drawn ROI must be within 10 pixels of the true region boundaries");

	}

	int[] wList = WindowManager.getIDList();
	String[] titles = new String[wList.length + 1];
	titles[0] = "No covariates";
	for (int i = 0; i < wList.length; i++) {
	    ImagePlus imp = WindowManager.getImage(wList[i]);
	    if (imp != null)
		titles[i + 1] = (imp).getTitle();
	    else
		titles[i + 1] = "";
	}
	gd.addChoice("Covariate image", titles, titles[0]);

	return gd;

    }

    static float[][] deepCopy(float[][] original) {
	float[][] copy = new float[original.length][];
	for (int i = 0; i < original.length; i++)
	    copy[i] = original[i].clone();
	return copy;
    }

}
