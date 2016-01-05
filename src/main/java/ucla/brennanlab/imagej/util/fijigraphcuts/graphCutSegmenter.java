package ucla.brennanlab.imagej.util.fijigraphcuts;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ucla.brennanlab.imagej.util.levelsets.GenericLevelSet;
import ucla.brennanlab.imagej.util.levelsets.RoiImageLevelSet;
import ucla.brennanlab.imagej.util.levelsets.shapeKernelDensityEstimation;
import ucla.brennanlab.imagej.util.stochastic.ImageStats;

/**
 * This class segments a given ImageProcessor. We will interpret the foreground
 * as source, and the background as the sink
 * 
 * @author Josh Chang joshchang@ucla.edu
 * 
 */
public class graphCutSegmenter {
    int width;
    int height;
    ImageProcessor ip;
    public GraphCut gc;
    public double lengthPenalty = 20;
    public float Energy;
    private ImageStats s;
    private shapeKernelDensityEstimation shapeKernelDensityEstimate;

    public boolean[][] innermask;
    public boolean[][] outermask;

    public graphCutSegmenter(ImageProcessor ip, double MU) {
	this(ip);
	this.lengthPenalty = MU;

    }

    public void setInnerMask(boolean[][] im) {
	this.innermask = im;
    }

    public graphCutSegmenter(ImageProcessor ip) {
	this.ip = ip;
	this.width = ip.getWidth();
	this.height = ip.getHeight();
	/**
	 * Initialize eight-neighbor graph cut neighborhood should be
	 * 4*(width-1)*(height-1)+width+height-2
	 */
	this.gc = new GraphCut(width * height, 4 * (width - 1) * (height - 1)
		+ width + height - 1);
    }

    /**
     * Return an Roi corresponding to the current graph cut This function uses
     * the LevelSet class for laziness, should write a method for this guy but
     * whatever
     * 
     * @return Roi corresponding to the boundary
     */
    // TODO: Don't be so lazy
    public Roi returnRoi() {
	boolean[][] mask = new boolean[width][height];
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		mask[x][y] = gc.getTerminal(y * width + x) == Terminal.FOREGROUND ? true
			: false;
	    }
	}
	RoiImageLevelSet rp = new RoiImageLevelSet(mask);

	return rp.getRoi(false);

    }

    public float fractionInner() {
	int inner = 0;
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		if (gc.getTerminal(y * width + x) == Terminal.FOREGROUND)
		    inner++;
	    }
	}
	return (float) (1.0 * inner / width / height);
    }

    public boolean[][] returnMask() {
	boolean[][] mask = new boolean[width][height];
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		mask[x][y] = gc.getTerminal(y * width + x) == Terminal.FOREGROUND ? true
			: false;
	    }
	}
	return mask;
    }

    public GenericLevelSet getLevelSet() {
	return new GenericLevelSet(returnMask());
    }

    public ImageProcessor returnMaskProcessor() {
	int[][] mask = new int[width][height];
	ByteProcessor ip1 = new ByteProcessor(width, height, new byte[width
		* height], null);
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		mask[x][y] = gc.getTerminal(y * width + x) == Terminal.FOREGROUND ? 0
			: 255;
		ip1.putPixel(x, y, mask[x][y]);
	    }
	}
	return (ip1);
    }

    public ImagePlus returnMaskedImage() {

	return (new ImagePlus("", returnMaskProcessor()));
    }

    public void setLengthPenalty(double MU) {
	this.lengthPenalty = MU;
    }

    /**
     * Set edge weights with no prior information other than length penalty
     * could just set a pointer maybe
     */
    public void setEdgeWeights() {
	for (int y = 0; y < height - 1; y++) {
	    for (int x = 0; x < width - 1; x++) {
		gc.setEdgeWeight(y * width + x, y * width + x + 1,
			(float) (lengthPenalty * Math.PI / 8));
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x,
			(float) (lengthPenalty * Math.PI / 8));
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x + 1,
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)));
		gc.setEdgeWeight(y * width + x + 1, (y + 1) * width + x,
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)));

	    }
	    gc.setEdgeWeight(y * width + width - 1,
		    (y + 1) * width + width - 1, (float) (lengthPenalty
			    * Math.PI / 8));
	}
	for (int x = 0; x < width - 1; x++) {
	    gc.setEdgeWeight((height - 1) * width + x, (height - 1) * width + x
		    + 1, (float) (lengthPenalty * Math.PI / 8));
	}
    }

    /**
     * The prior shape is a reference shape, or rather a prior guess at the true
     * segmentation that is used for determining weighting according to the MM
     * algorithm. Setting these weights should be parallelizable
     * 
     * @param skde
     *            Shape kernel density estimation
     * @param priorShape
     *            The segmentation guess
     */
    public void setEdgeWeights(shapeKernelDensityEstimation skde,
	    GenericLevelSet priorShape) {
	// what I would like to do is reset the shape weights rather
	// than recompute and reset edge weights
	gc.resetEdgeNum();
	float weight;
	/**
	 * Calculate contributions from the density estimation
	 */

	double[] kernelWeights = new double[skde.shapes.size()];

	int j = 0;
	for (; j < skde.shapes.size(); j++) {
	    kernelWeights[j] = Math
		    .exp(-0.5
			    * Math.log(2 * Math.PI * skde.sigmaSquared)
			    - 0.5
			    * (skde.computeDistance(priorShape, skde.shapes
				    .get(j)) / skde.sigmaSquared))
		    * skde.weights[j];

	}
	kernelWeights = normalize(kernelWeights);
	GenericLevelSet currentshape;
	// weights are length penalty + shape penalty
	float[] weights = new float[] { 0, 0, 0, 0 };

	for (int y = 0; y < height - 1; y++) {
	    for (int x = 0; x < width - 1; x++) {
		weights = new float[] { (float) (lengthPenalty * Math.PI / 8),
			(float) (lengthPenalty * Math.PI / 8),
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)),
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)) };

		for (j = 0; j < skde.shapes.size(); j++) {
		    currentshape = skde.shapes.get(j);
		    weights[0] += 0.5
			    * kernelWeights[j]
			    * Math.pow(Math.abs(0.5 * currentshape.get(x, y)
				    + 0.5 * currentshape.get(x + 1, y)), skde
				    .getAlpha()) / skde.sigmaSquared;
		    weights[1] += 0.5
			    * kernelWeights[j]
			    * Math.pow(
				    Math.abs(0.5 * currentshape.get(x, y) + 0.5
					    * currentshape.get(x, y + 1))
					    * Math.abs(0.5
						    * currentshape.get(x, y)
						    + 0.5
						    * currentshape
							    .get(x, y + 1)),
				    skde.getAlpha()) / skde.sigmaSquared;
		    weights[2] += kernelWeights[j]
			    * Math.pow(
				    Math.abs(0.5 * currentshape.get(x, y) + 0.5
					    * currentshape.get(x + 1, y + 1))
					    * Math.abs(0.5
						    * currentshape.get(x, y)
						    + 0.5
						    * currentshape
							    .get(x, y + 1)),
				    skde.getAlpha()) / skde.sigmaSquared
			    / Math.sqrt(2);
		    weights[3] += kernelWeights[j]
			    * Math.pow(
				    Math.abs(0.5 * currentshape.get(x + 1, y)
					    + 0.5 * currentshape.get(x, y + 1))
					    * Math.abs(0.5
						    * currentshape.get(x, y)
						    + 0.5
						    * currentshape
							    .get(x, y + 1)),
				    skde.getAlpha()) / skde.sigmaSquared
			    / Math.sqrt(2);
		}
		gc.setEdgeWeight(y * width + x, y * width + x + 1, weights[0]);
		gc
			.setEdgeWeight(y * width + x, (y + 1) * width + x,
				weights[1]);
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x + 1,
			weights[2]);
		gc.setEdgeWeight(y * width + x + 1, (y + 1) * width + x,
			weights[3]);

	    }

	    weight = (float) (lengthPenalty * Math.PI / 8);
	    for (j = 0; j < skde.shapes.size(); j++) {

		weight += 0.5
			* kernelWeights[j]
			* Math.pow(Math.abs(0.5
				* skde.shapes.get(j).get(width - 1, y) + 0.5
				* skde.shapes.get(j).get(width - 1, y + 1)),
				skde.getAlpha()) / skde.sigmaSquared;

	    }
	    gc.setEdgeWeight(y * width + width - 1,
		    (y + 1) * width + width - 1, weight);
	}
	for (int x = 0; x < width - 1; x++) {
	    weight = (float) (lengthPenalty * Math.PI / 8);
	    for (j = 0; j < skde.shapes.size(); j++) {

		weight += 0.5
			* kernelWeights[j]
			* Math.pow(Math.abs(0.5
				* skde.shapes.get(j).get(x, height - 1) + 0.5
				* skde.shapes.get(j).get(x + 1, height - 1)),
				skde.getAlpha()) / skde.sigmaSquared;

	    }
	    gc.setEdgeWeight((height - 1) * width + x, (height - 1) * width + x
		    + 1, weight);
	}
    }

    /**
     * Set the graph cut weights
     * 
     * @param prior
     *            The prior probability of being in foreground
     * @param s
     *            The image statistics
     * @param skde
     *            Shape kernel density estimate
     */

    public void setNodeWeights(shapeKernelDensityEstimation skde,
	    GenericLevelSet priorShape, ImageStats s) {

	double[] kernelWeights = new double[skde.shapes.size()];

	for (int j = 0; j < skde.shapes.size(); j++) {
	    kernelWeights[j] = skde.weights[j]
		    * Math.exp(-0.5
			    * Math.log(2 * Math.PI * skde.sigmaSquared)
			    - 0.5
			    * (skde.computeDistance(priorShape, skde.shapes
				    .get(j)) / skde.sigmaSquared));

	    if (Double.isNaN(kernelWeights[j])) {
		kernelWeights[j] = 0;
	    }
	}

	final double[] weights = normalize(kernelWeights);

	float source, sink;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {

		source = 0;
		sink = 0;

		double pixelval = this.ipGetSafeValue(x, y);
		// just ignore outliers
		if (Math.abs(pixelval) > s.innerinnermean + 3 * s.innerinnersd) {

		    pixelval = (s.innerinnermean / s.innerinnersd + s.outermean
			    / s.outersd)
			    / (1.0 / s.innerinnersd + 1.0 / s.outersd);
		}
		// likelihood portion of the weight

		if (Double.isNaN(source) || Double.isNaN(sink))
		    IJ.log("Kernel object sigma^2: " + skde.sigmaSquared);
		
		// Hack for blood vessels, if a pixel is 0, and 3 of its neighbors are zero
		// don't add it to the likelihood
		
		int deadNeighbors = 0;
		if(pixelval==0){
		    if(this.ipGetSafeValue(x, y+1)==0) deadNeighbors++;
		    if(this.ipGetSafeValue(x, y-1)==0) deadNeighbors++;
		    if( this.ipGetSafeValue(x+1, y)==0) deadNeighbors++;
		    if( this.ipGetSafeValue(x-1, y)==0) deadNeighbors++;
		    if( this.ipGetSafeValue(x-1, y-1)==0) deadNeighbors++;
		    if( this.ipGetSafeValue(x+1, y-1)==0) deadNeighbors++;
		    if( this.ipGetSafeValue(x-1, y+1)==0) deadNeighbors++;
		    if( this.ipGetSafeValue(x+1, y+1)==0) deadNeighbors++;
		}
		
		if(deadNeighbors<4) {
		    source += 0.5
			    * Math.log(2 * Math.PI * s.outersd * s.outersd)
			    + 0.5
			    * Math.pow((pixelval - s.outermean) / s.outersd, 2);
		    sink += 0.5
			    * Math.log(2 * Math.PI * s.innerinnersd
				    * s.innerinnersd)
			    + 0.5
			    * Math.pow((pixelval - s.innerinnermean)
				    / s.innerinnersd, 2);
		}
		for (int j = 0; j < skde.shapes.size(); j++) {
		    if (skde.shapes.get(j).get(x, y) < 0.0) {
			source += weights[j]
				* Math.pow(Math.abs(skde.shapes.get(j)
					.get(x, y)), skde.getAlpha()) / 2
				/ skde.sigmaSquared;
		    } else
			sink += weights[j]
				* Math.pow(Math.abs(skde.shapes.get(j)
					.get(x, y)), skde.getAlpha()) / 2
				/ skde.sigmaSquared;
		}

		// Near image edges, adjust

		if (x == 0 || y == 0 || x == width - 1 || y == height - 1) {
		    sink += this.lengthPenalty * Math.PI
			    * (1 + 2 * Math.pow(2, -0.5)) / 8;
		    for (int j = 0; j < skde.shapes.size(); j++) {
			sink += weights[j]
				* Math.pow(Math.abs(skde.shapes.get(j)
					.get(x, y) + 0.5), skde.getAlpha()) / 2
				/ skde.sigmaSquared;

		    }
		}

		if (Double.isNaN(source) || Double.isNaN(sink)) {
		    source = 0;
		    sink = 0;
		}
		gc.setTerminalWeights(y * width + x, source, sink);
	    }
	}

    }
    
    public double ipGetSafeValue(int x, int y){
	x = x<0 ? 0 : x;
	x = x>=width ? width-1 : x;
	y = y<0 ? 0 : y;
	y= y>=height ? height-1: y;
	return ip.getPixelValue(x, y);
	
    }

    /**
     * Chan - Vese node weights
     */
    public void setNodeWeights(ImageStats s) {
	this.setS(s);

	float source, sink;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		source = 0;
		sink = 0;

		double pixelval = ip.getPixelValue(x, y);
		if (Double.isNaN(source) || Double.isNaN(sink))
		    IJ.showMessage("BLAH1");
		source += 0.5 * Math.log(2 * Math.PI * s.outersd * s.outersd)
			+ 0.5
			* Math.pow((pixelval - s.outermean) / s.outersd, 2);
		sink += 0.5 * Math.log(2 * Math.PI * s.innersd * s.innersd)
			+ 0.5
			* Math.pow((pixelval - s.innermean) / s.innersd, 2);

		if (Double.isNaN(source) || Double.isNaN(sink)) {
		    source = 0;
		    sink = 0;

		}

		gc.setTerminalWeights(y * width + x, source, sink);
	    }
	}

    }

    public float relaxEnergy() {
	this.Energy = gc.computeMaximumFlow(true, null);
	return (this.Energy);
    }

    public float overlap(boolean[][] mask1, boolean[][] mask2) {
	int w1 = mask1.length;
	int w2 = mask2.length;
	int h1 = mask1[0].length;
	int h2 = mask2[0].length;
	if (w1 != w2 || h1 != h2)
	    return 0;
	float area = 0;
	int over = 0;
	int in = 0;
	for (int x = 0; x < w1; x++) {
	    for (int y = 0; y < h1; y++) {
		if (!mask1[x][y]) {
		    in++;
		    if (!mask2[x][y]) {
			over++;
		    }

		}
	    }
	}
	area = (float) over / in;
	return area;

    }

    public double[] normalize(double[] numbers) {
	if (numbers.length == 1)
	    return new double[] { 1 };
	double sum = 0;
	double[] num = new double[numbers.length];
	for (int i = 0; i < numbers.length; i++) {
	    sum += numbers[i];
	}
	for (int i = 0; i < numbers.length; i++) {
	    num[i] = sum == 0 ? (double) 1 / numbers.length : numbers[i] / sum;
	}
	return num;
    }

    public void setS(ImageStats s) {
	this.s = s;
    }

    public ImageStats getS() {
	return s;
    }

    public void setShapeKernelDensityEstimate(
	    shapeKernelDensityEstimation shapeKernelDensityEstimate) {
	this.shapeKernelDensityEstimate = shapeKernelDensityEstimate;
    }

    public shapeKernelDensityEstimation getShapeKernelDensityEstimate() {
	return shapeKernelDensityEstimate;
    }

}
