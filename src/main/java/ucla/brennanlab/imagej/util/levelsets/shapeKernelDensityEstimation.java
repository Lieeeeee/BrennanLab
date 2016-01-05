package ucla.brennanlab.imagej.util.levelsets;

import java.util.ArrayList;
import java.util.Iterator;

public class shapeKernelDensityEstimation {
    public ArrayList<GenericLevelSet> shapes;
    public double[] weights;
    public double sigmaSquared;
    public double alpha = 1;
    int width, height;

    public boolean[][] strictlyInBackground;
    public boolean[][] strictlyInForeground;

    public float[][] pixelProbability; // 
    public boolean[][] possibleForegroundMask;

    /**
     * Hard constraint
     * 
     * @param inadmissibleRegion
     */
    public void setForeGround(boolean[][] strictlyInForeground) {
	this.strictlyInForeground = strictlyInForeground;
    }

    /**
     * Hard constraint
     * 
     * @param inadmissibleRegion
     */
    public void setBackGround(boolean[][] strictlyInBackground) {
	this.strictlyInBackground = strictlyInBackground;
    }

    /**
     * Take the shapes, sort them, compute \sigma^2
     * 
     * @param weights
     * @param shapes
     */
    public shapeKernelDensityEstimation(double[] weights,
	    ArrayList<GenericLevelSet> shapes) {
	this.shapes = shapes;
	this.weights = weights;
	this.width = this.shapes.get(0).width;
	this.height = this.shapes.get(0).height;
	if (shapes.size() != weights.length)
	    return;
	/**
	 * compute the pixel probabilities, for regularization perhaps
	 */
	this.pixelProbability = new float[width][height];
	
    }

    public void estimateSigmaSquared() {

	/**
	 * Go through shapes, calculate distances
	 */
	double minDistances[] = new double[this.shapes.size()];
	double tempDistance;
	double sigmaS = 0;
	for (int j = 0; j < this.shapes.size(); j++) {
	    minDistances[j] = Double.MAX_VALUE;
	    for (int k = 0; k < this.shapes.size(); k++) {
		if (k == j)
		    continue;
		else {
		    tempDistance = computeDistance(shapes.get(k), shapes.get(j));
		    if (tempDistance < minDistances[j])
			minDistances[j] = tempDistance;
		}
	    }
	    sigmaS += this.weights[j] * minDistances[j];
	}
	this.sigmaSquared = sigmaS;
    }

    /**
     * Compute the non-symmetric shape distance. Penalize edges in shape 1, by
     * the corresponding signed distance from shape 2.
     * 
     * 
     * @param shape1
     *            Shape to compare
     * @param shape2
     *            Acts as reference shape
     * @return The distance between the shapes
     */
    public static double computeL1Distance(final GenericLevelSet shape1,
	    final GenericLevelSet shape2) {
	double d = 0;
	/**
	 * First compute the pixel-wise differences
	 */

	for (int x = 0; x < shape1.width; x++) {
	    for (int y = 0; y < shape1.height; y++) {
		if (shape1.get(x, y) * shape2.get(x, y) < 0) {
		    d += Math.pow(Math.abs(shape2.get(x, y)), 1);
		}
	    }
	}

	/**
	 * Now compute the pair-wise differences
	 */
	for (int y = 0; y < shape1.height - 1; y++) {
	    for (int x = 0; x < shape1.width - 1; x++) {
		if (shape1.get(x, y) * shape1.get(x + 1, y) < 0)
		    d += Math.pow(Math.abs(0.5 * shape2.get(x, y) + 0.5
			    * shape2.get(x + 1, y)), 1);
		if (shape1.get(x, y) * shape1.get(x, y + 1) < 0)
		    d += Math.pow(Math.abs(0.5 * shape2.get(x, y) + 0.5
			    * shape2.get(x, y + 1)), 1);
		if (shape1.get(x + 1, y + 1) * shape1.get(x, y) < 0)
		    d += Math.sqrt(2)
			    * Math.pow(Math.abs(0.5 * shape2.get(x, y) + 0.5
				    * shape2.get(x + 1, y + 1)), 1);
		if (shape1.get(x + 1, y) * shape1.get(x, y + 1) < 0)
		    d += Math.sqrt(2)
			    * Math.pow(Math.abs(0.5 * shape2.get(x + 1, y)
				    + 0.5 * shape2.get(x, y + 1)), 1);
	    }
	}
	return d;

    }

    /**
     * Compute the non-symmetric shape distance. Penalize edges in shape 1, by
     * the corresponding signed distance from shape 2.
     * 
     * 
     * @param shape1
     *            Shape to compare
     * @param shape2
     *            Acts as reference shape
     * @return
     */
    public double computeDistance(GenericLevelSet shape1, GenericLevelSet shape2) {
	double d = 0;
	/**
	 * First compute the pixel-wise differences
	 */

	for (int x = 0; x < shape1.width; x++) {
	    for (int y = 0; y < shape1.height; y++) {
		if (shape1.get(x, y) * shape2.get(x, y) < 0) {
		    d +=  Math.pow(Math.abs(shape2.get(x, y)), alpha);
		}
	    }
	}

	/**
	 * Now compute the pair-wise differences
	 */
	for (int y = 0; y < shape1.height - 1; y++) {
	    for (int x = 0; x < shape1.width - 1; x++) {
		if (shape1.get(x, y) * shape1.get(x + 1, y) < 0)
		    d += Math.pow(Math.abs(0.5 * shape2.get(x, y) + 0.5
			    * shape2.get(x + 1, y)), alpha);
		if (shape1.get(x, y) * shape1.get(x, y + 1) < 0)
		    d += Math.pow(Math.abs(0.5 * shape2.get(x, y) + 0.5
			    * shape2.get(x, y + 1)), alpha);
		if (shape1.get(x + 1, y + 1) * shape1.get(x, y) < 0)
		    d += Math.sqrt(2)
			    * Math.pow(Math.abs(0.5 * shape2.get(x, y) + 0.5
				    * shape2.get(x + 1, y + 1)), alpha);
		if (shape1.get(x + 1, y) * shape1.get(x, y + 1) < 0)
		    d += Math.sqrt(2)
			    * Math.pow(Math.abs(0.5 * shape2.get(x + 1, y)
				    + 0.5 * shape2.get(x, y + 1)), alpha);
	    }
	}
	return d;

    }

    /**
     * Use the KDE to compute the probability density
     * 
     * @param shape
     * @return
     */
    public double estimateDensity(GenericLevelSet shape) {
	double prob = 0;
	double[] normalizedDistances = new double[shapes.size()];
	Iterator<GenericLevelSet> it = shapes.iterator();
	int i = 0;
	while (it.hasNext()) {
	    normalizedDistances[i] = computeDistance(shape, it.next())
		    / sigmaSquared;
	    prob += Math.exp(-0.5 * Math.log(2 * Math.PI * sigmaSquared) - 0.5
		    * normalizedDistances[i]);
	    i++;
	}
	return prob;
    }

    public void setAlpha(double alpha) {
	this.alpha = alpha;
    }

    public double getAlpha() {
	return alpha;
    }
    
   
}
