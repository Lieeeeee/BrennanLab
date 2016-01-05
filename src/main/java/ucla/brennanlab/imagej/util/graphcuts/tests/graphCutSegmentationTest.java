package ucla.brennanlab.imagej.util.graphcuts;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import ucla.brennanlab.imagej.util.fijigraphcuts.graphCutSegmenter;
import ucla.brennanlab.imagej.util.stochastic.ImageStats;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class graphCutSegmentationTest implements PlugInFilter {

    ImagePlus imp;

    int repetitions = 10;

    public void run(ImageProcessor ip) {
	NumberFormat fmt = new DecimalFormat("####.##");
	IJ
		.log("Welcome to the Naive Chan-Vese-like graph cuts segmentation test");
	graphCutSegmenter gcs = new graphCutSegmenter(ip);
	double lengthpenalty=10;
	
	GenericDialog gd = new GenericDialog("Choose the length penalty");
	gd.addNumericField("Length penalty", lengthpenalty, 2);
	gd.showDialog();
	if (gd.wasCanceled()) {
	    return;
	}
	
	gcs.setLengthPenalty(gd.getNextNumber());
	IJ.log("Computing initial image stats");
	ImageStats s = new ImageStats();

	boolean[][] mask = new boolean[ip.getWidth()][ip.getHeight()];
	ImageProcessor ip2 = ip.duplicate();
	GaussianBlur gb = new GaussianBlur();
	gb.blur(ip2, 2);
	ImageStatistics is = ImageStatistics.getStatistics(ip2,
		ij.measure.Measurements.MEAN + ij.measure.Measurements.STD_DEV,
		null);

	double mean = is.mean;
	double sd = is.stdDev;

	for (int x = 0; x < ip.getWidth(); x++) {
	    for (int y = 0; y < ip.getHeight(); y++) {
		mask[x][y] = ip2.getPixelValue(x, y) > mean + sd ? true : false;
	    }
	}

	s.updateStats(ip2, mask);
	IJ.log("Done computing initial image stats");
	RoiManager rm = RoiManager.getInstance();
	if (rm != null)
	    rm.setVisible(true);
	else
	    rm = new RoiManager();

	ImagePlus imp1 = imp.createImagePlus();
	imp1.setProcessor("", (FloatProcessor) ip.duplicate());
	imp1.show();
	gcs.setEdgeWeights();
	int i = 1;

	while (i < repetitions) {
	    IJ.log("Mean inner intensity: " + fmt.format(s.innermean)
		    + "  outer intensity: " + fmt.format(s.outermean));
	    IJ.log("Starting graph-cuts segmentation iteration " + i);
	    long startTime = System.nanoTime();
	    long endTime;
	    gcs.gc.reset();
	    gcs.setNodeWeights(s);
	    
	    try {
		IJ.log("Energy "+gcs.relaxEnergy());
	    } finally {
		endTime = System.nanoTime();

	    }
	     Roi roi =gcs.returnRoi();
	     if(roi==null) IJ.showMessage("woof");
	    s.updateStats(ip2, gcs.returnMask());
	    roi.setName("Iteration " + i);
	    roi.setStrokeWidth(2);
	    rm.addRoi((Roi) roi.clone());
	
	    IJ.log("Finished repetition " + i + " in "
		    + fmt.format((double) (endTime - startTime) / 1000000)
		    + " ms");
	   
	    i++;
	}

	Roi r = gcs.returnRoi();
	r.setStrokeWidth(3);
	imp1.getProcessor().draw(r);
	imp1.setTitle("Segmentation results");
	imp1.updateAndDraw();

    }

    public int setup(String arg, ImagePlus imp) {
	this.imp = imp;
	return DOES_8G + DOES_16 + DOES_32 + NO_CHANGES + CONVERT_TO_FLOAT;
    }

}