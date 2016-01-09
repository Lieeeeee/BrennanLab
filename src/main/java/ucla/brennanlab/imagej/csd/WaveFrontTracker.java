package ucla.brennanlab.imagej.csd;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;
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
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ucla.brennanlab.imagej.graphcuts.GraphCutSegmenter;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;
import ucla.brennanlab.imagej.util.levelsets.ShapePrior;
import ucla.brennanlab.imagej.util.stochastic.GaussianMixture;
import ucla.brennanlab.imagej.util.stochastic.IntensityModel;
import ucla.brennanlab.imagej.util.stochastic.LaplaceMixture;
import ucla.brennanlab.imagej.util.stochastic.SpeedField;

import java.awt.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;


/**
 * Boundary segmentation using the kriging speed model
 *
 * @author Josh Chang
 */
public class WaveFrontTracker implements PlugInFilter {
    /**
     * the number of importance samples is J*(2L+1). L is usually the dimension
     * of the stochastic state, here it will be the number of points at which we
     * will estimate the speed, interpolating speeds using these points
     */
    ImagePlus imp;
    int timesteps;
    int slices;
    int maxiters = 10;
    double initialLengthPenalty = 1;
    double lengthPenalty = 1;
    double alpha = 2;
    boolean showPredictions = false;
    boolean useManualInitializationPrior = false;
    Roi userDrawnRoi;
    int speedSamples = 16;
    double priorSpeedMean = 5, priorSpeedSD = 3;
    int width;
    int height;
    int innerBandWidth = 25;
    SpeedField maxAposteriorSpd; // filtered speed field
    ArrayList<Roi> roisOfSegmentations;
    boolean csdHasStarted = false;
    Color roiColor = new Color(255, 105, 180, 255); // Color to draw the ROI
    final String[] likelihoodmodels = {"Gaussian", "Laplacian"};
    String likelihoodmodel;
    String covariateImage;

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

        double inputspeed = gd.getNextNumber();
        double inputstder = gd.getNextNumber();
        alpha = gd.getNextNumber();
        lengthPenalty = gd.getNextNumber();
        initialLengthPenalty = gd.getNextNumber();
        priorSpeedMean = inputspeed;
        priorSpeedSD = inputstder;
        this.speedSamples = (int) gd.getNextNumber();
        this.innerBandWidth = (int) gd.getNextNumber();
        this.showPredictions = gd.getNextBoolean();
        this.likelihoodmodel = gd.getNextChoice();
        this.covariateImage = gd.getNextChoice();


        if (userDrawnRoi != null)
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
        RoiManager roiman = RoiManager.getInstance();
        if (roiman != null)
            roiman.setVisible(true);
        else
            roiman = new RoiManager();

        GraphCutSegmenter gcSegmenter = new GraphCutSegmenter(ip);

        ImageStack gSegmentStack = imp.createEmptyStack();
        ImagePlus gSegmentImp = imp.createImagePlus();
        gSegmentImp.setTitle(imp.getTitle() + "_segmentation");

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
                this.lengthPenalty = imp.getWidth()/40.0;
            }

            if (csdHasStarted) {

                ImplicitShape2D prevLS = new ImplicitShape2D(gcSegmenter.returnMask());
                if (gcSegmenter.fractionInner() < 10.0 / width / height) {
                    IJ.log("Less than 10 pixels in the inner region");
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
                meanNext.advectUniformSpeed(this.priorSpeedMean, 1.0);

                // Infer the speed field here!!! @TODO

                Roi nextRoi = meanNext.getRoi(false);
                imp.setRoi(nextRoi);

                nextRoi.setName("Predicted mean position " + currentSlice);
                roiman.addRoi(nextRoi);

                double[] effectiveWeights = new double[this.speedSamples];
                double speed;
                for (int i = 0; i < this.speedSamples; i++) {
                    speed = normalRV.nextDouble(priorSpeedMean, priorSpeedSD);

                    positions.add(new ImplicitShape2D(
                            prevLS.getMask()));
                    effectiveWeights[i] = 1;// redo the weighting soon

                    positions.get(i).advectUniformSpeed(speed, 1.0);

                    Roi pred = positions.get(i).getRoi(true);

                    if (pred == null)
                        continue;

                    if (showPredictions) {
                        pred
                                .setName("Sample " + (i + 1)
                                        + " of predicted position of "
                                        + (currentSlice));
                        roiman.addRoi(pred);
                    }

                }

                effectiveWeights = normalize(effectiveWeights);

                /*************************************************************
                 * Kernel density estimation on the shapes to determine prior
                 *************************************************************/

                ShapePrior shapePriorDensity = new ShapePrior(
                        positions, effectiveWeights);

                IJ.log("Prediction distance variance " + shapePriorDensity.beta);

                IJ.log("... done computing priors");
                IJ.log("... computing maxflow for slice " + currentSlice
                        + "...");
                long startTime = System.nanoTime();
                float minE;
                float tempE;
                long endTime;

                IntensityModel likelihood;
                if(this.likelihoodmodel=="Gaussian") {
                    likelihood = new GaussianMixture();
                }else{
                    likelihood = new LaplaceMixture();
                }

                likelihood.Infer(currentImageProcessor, meanNext
                        .getMask());

                gcSegmenter = new GraphCutSegmenter(currentImageProcessor);
                gcSegmenter.setLengthPenalty((float) this.lengthPenalty);
                gcSegmenter.setIntensityModel(likelihood);

                /********************************************************************************
                 * Segment NOW
                 *********************************************************************************/

                float currentMeanSpeed =  (float)this.maxAposteriorSpd.latestSpeed;

                try {
                    gcSegmenter.setNodeWeights(shapePriorDensity, meanNext, likelihood);
                    gcSegmenter.setEdgeWeights(shapePriorDensity, meanNext);
                    minE = gcSegmenter.relaxEnergy();
                    /**
                     * MM iterations
                     */
                    likelihood.Infer(currentImageProcessor, gcSegmenter
                            .returnMask());
                    gcSegmenter.setEdgeWeights(shapePriorDensity, gcSegmenter.getLevelSet());
                    long currentTime;
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
                        gcSegmenter.gc.reset(); //  Don't remember what this does

                        gcSegmenter.setNodeWeights(shapePriorDensity, gcSegmenter.getLevelSet(), likelihood);
                        gcSegmenter.setEdgeWeights(shapePriorDensity, gcSegmenter.getLevelSet());
                        tempE = gcSegmenter.relaxEnergy();
                        currentTime = System.nanoTime();
                        IJ.log("\t     Iter " + i + " Elapsed time: "
                                + ((float) (currentTime - startTime) / 1000000)
                                + " ms" + " Energy " + tempE);
                        if (minE - tempE < 1 && i > 2 || Double.isNaN(minE)
                                || Double.isInfinite(minE) || i> this.maxiters) {
                            minE = tempE;
                            if (Double.isNaN(minE) || Double.isInfinite(minE)) {
                                IJ.log("Encountered invalid energy");
                            }

                            break;
                        }
                        minE = tempE;
                        // if (i % 3 == 0) {
                        try {
                            Roi currentProgress = (Roi) gcSegmenter.getLevelSet().getRoi(
                                    false).clone();
                            currentProgress.setStrokeColor(roiColor);

                            currentProgress.setName("Position_" + currentSlice
                                    + "_iteration_" + i);
                            imp.setRoi(currentProgress);
                            roiman.addRoi(currentProgress);

                        } catch (NullPointerException e) {
                            csdHasStarted = false;

                        }

                        // }
                        i++;
                    }

                } finally {
                    endTime = System.nanoTime();
                }
                Roi segmentedRoi = gcSegmenter.returnRoi();
                segmentedRoi.setStrokeColor(roiColor);


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

                    gSegmentImp.setStack(gSegmentStack);
                    gSegmentImp.show();
                    gSegmentImp.setSlice(gSegmentImp.getStackSize());
                    gSegmentImp.updateAndDraw();

                    roisOfSegmentations.add((Roi) segmentedRoi.clone());
                    roisOfSegmentations.get(currentSlice - firstCSDslice)
                            .setName("Position_" + currentSlice);
                    roisOfSegmentations.get(currentSlice - firstCSDslice)
                            .setStrokeColor(this.roiColor);
                    roisOfSegmentations.get(currentSlice - firstCSDslice)
                            .setStrokeWidth(3);

                    roiman.addRoi(roisOfSegmentations.get(currentSlice
                            - firstCSDslice));


                    gSegmentStack.getProcessor(currentSlice - firstCSDslice + 1).setColor(this.roiColor);
                    (gSegmentStack
                            .getProcessor(currentSlice - firstCSDslice + 1))
                            .drawRoi(
                                    roisOfSegmentations.get(currentSlice
                                            - firstCSDslice));
                    gSegmentImp.setSlice(currentSlice - firstCSDslice + 1);

                } else {
                    IJ.log("Something went wrong, I was unable to obtain a segmentation...");

                    lastCSDslice = currentSlice - 1;
                    continue;
                }

                long duration = endTime - startTime;
                double dist = shapePriorDensity.computeDistance(gcSegmenter.getLevelSet(), meanNext)
                        * shapePriorDensity.beta;
                IJ.log("Computed maxflow energy of " + minE + " in "
                        + (float) duration / 1000000 + " ms");
                IJ
                        .log("Normalized distance between mean prediction and result: "
                                + dist);

                if (Double.isNaN(minE) || minE == 0 || Double.isNaN(dist)) {
                    lastCSDslice = currentSlice - 1;
                    continue;
                }

                ImplicitShape2D currentPosition = gcSegmenter.getLevelSet();
                this.maxAposteriorSpd.addArrival(currentPosition);
                this.maxAposteriorSpd.updateGMRF();

                currentMeanSpeed = (float) this.maxAposteriorSpd.latestSpeed;

                IJ.log("Current average wavespeed is "
                        + currentMeanSpeed + " pixels per frame"); // current wavespeed

            } else {
                /****************************************************
                 * CSD maybe hasn't started yet
                 *****************************************************/

                if (this.useManualInitializationPrior) {
                    this.userDrawnRoi.setName("User_drawn_template_slice_"
                            + currentSlice);
                    roiman.addRoi(this.userDrawnRoi);
                    ImplicitShape2D manualLS = (new ImplicitShape2D(
                            this.userDrawnRoi, width, height));
                    ArrayList<ImplicitShape2D> userDefinedShapes = new ArrayList<ImplicitShape2D>(
                            1);
                    userDefinedShapes.add(manualLS);
                    double[] weights = new double[1];
                    weights[0] = 1;
                    ShapePrior userKDE = new ShapePrior(
                            userDefinedShapes, weights);
                    userKDE.setMultiplier(1.0);
                    userKDE.beta = 1.0;

                    IntensityModel s;
                    if(this.likelihoodmodel=="Gaussian") {
                        s = new GaussianMixture();
                    }else{
                        s = new LaplaceMixture();
                    }
                    s.Infer(currentImageProcessor, manualLS.getMask());

                    gcSegmenter = new GraphCutSegmenter(currentImageProcessor, 1);

                    gcSegmenter.setIntensityModel(s);
                    gcSegmenter.setEdgeWeights(userKDE, manualLS);
                    gcSegmenter.setNodeWeights(userKDE, manualLS, s);
                    double minE = gcSegmenter.relaxEnergy();
                    IJ.log(minE + " " + gcSegmenter.fractionInner());
                    s.Infer(currentImageProcessor, gcSegmenter.returnMask());

                } else {
                    IJ.log("... computing naive maxflow for slice "
                            + currentSlice + " ... ");
                    long startTime = System.nanoTime();
                    long endTime;
                    try {
                        gcSegmenter = this.naivelyFindSegmentation(
                                currentImageProcessor, 3);

                    } finally {
                        endTime = System.nanoTime();
                    }

                    long duration = endTime - startTime;

                    IJ.log("Computed maxflow energy of " + gcSegmenter.Energy + " in "
                            + (float) duration / 1000000 + " ms");
                }
                // double inner = gcs.getS().innerarea;
                if (gcSegmenter.fractionInner() > 0.01) {
                    /**
                     * Provisionally say that the CSD has been initialized, but
                     * we might get kicked back here again if the subsequent
                     * frames yield inadmissible wave speed profiles
                     */
                    csdHasStarted = true;
                    IJ.log("Found a " + fmt.format(gcSegmenter.fractionInner() * 100)
                            + "% blob that might be the start of CSD at slice "
                            + currentSlice + "...");
                    Roi initialROI = gcSegmenter.returnRoi();
                    initialROI.setName("Position_" + currentSlice);
                    roiman.addRoi(initialROI);
                    firstCSDslice = currentSlice;
                    roisOfSegmentations.add((Roi) initialROI.clone());
                    this.maxAposteriorSpd.addArrival(gcSegmenter.getLevelSet());

                } else {
                    csdHasStarted = false;
                    IJ.log("Found no trace of CSD, skipping to next slice\n");
                }

            }

            currentSlice++;

        }

        IJ.setSlice(currentSlice);

        gSegmentImp.setOverlay(null);
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

        ImageProcessor speedIP = new FloatProcessor(speeds);
        speedIP.setMinAndMax(this.priorSpeedMean - this.priorSpeedSD / 2, this.priorSpeedMean * 1.5);
        ImagePlus speedImage = new ImagePlus("Speed_field_pixels_per_slice", speedIP);
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


        boolean[][] mask = new boolean[ip.getWidth()][ip.getHeight()];
        ImageProcessor ipCopy = ip.duplicate();
        GaussianBlur gb = new GaussianBlur();
        gb.blurGaussian(ipCopy, 3,3,0.01);

        ImageStatistics is = ImageStatistics.getStatistics(ipCopy,
                ij.measure.Measurements.MEAN + ij.measure.Measurements.STD_DEV,
                null);

        double mean = is.mean;
        double sd = is.stdDev;

        // Just some initial values
        for (int x = 0; x < ip.getWidth(); x++) {
            for (int y = 0; y < ip.getHeight(); y++) {
                mask[x][y] = ipCopy.getPixelValue(x, y) > mean + 2*sd;
            }
        }
        IntensityModel likelihood;
        if(this.likelihoodmodel=="Gaussian") {
            likelihood = new GaussianMixture();
        }else{
            likelihood = new LaplaceMixture();
        }
        likelihood.Infer(ipCopy, mask);

        GraphCutSegmenter gcSegmenter = new GraphCutSegmenter(ipCopy);
        gcSegmenter.setIntensityModel(likelihood);
        gcSegmenter.setLengthPenalty((float) this.initialLengthPenalty);
        int i = 1;
        gcSegmenter.setEdgeWeights((float)this.initialLengthPenalty);


        while (i++ < iterations) {
            IJ.log(""+i);
            gcSegmenter.setNodeWeights(likelihood);
            IJ.log("Energy " + gcSegmenter.relaxEnergy());
            likelihood.Infer(ipCopy, gcSegmenter.returnMask());

            IJ.log("\\mu_in " + likelihood.getPosteriorMean(true) + " \\mu_out"+ likelihood.getPosteriorMean(false));
            IJ.log("\\s_in " + likelihood.getPosteriorPrecision(true) + " \\s_out" + likelihood.getPosteriorPrecision(false));
        }
        return gcSegmenter;

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

        gd.addNumericField("Prior mean wave speed", this.priorSpeedMean, 2, 5,
                "pixels per frame");
        gd.addNumericField("standard error", this.priorSpeedSD, 2, 5,
                "pixels per frame");
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
        gd.addChoice("Likelihood model", likelihoodmodels, likelihoodmodels[0]);

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

    /**
     * Main method for debugging.
     *
     * For debugging, it is convenient to have a method that starts ImageJ, loads an
     * image and calls the plugin, e.g. after setting breakpoints.
     *
     * @param args unused
     */
    public static void main(String[] args) {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        Class<?> clazz = WaveFrontTracker.class;
        String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
        String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
        System.setProperty("plugins.dir", pluginsDir);

        // start ImageJ
        new IJ();

        // open the Clown sample
        ImagePlus image = IJ.openImage("http://imagej.net/images/clown.jpg");
        image.show();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }
}
