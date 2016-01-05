package ucla.brennanlab.imagej.util.levelsets;

import java.util.ArrayList;

import ucla.brennanlab.imagej.util.stochastic.SpeedField;

/**
 * Generic Level set class with to compute both the level set and Fast marching
 * methods
 * 
 * @author josh
 * 
 */
public class GenericLevelSet {

    public float[][] signedDistance;
    int width, height;
    static final double EPS = 1e-8;

    public GenericLevelSet(boolean[][] mask) {
	this.width = mask.length;
	this.height = mask[0].length;
	this.signedDistance = new float[width][height];
	int x = 0, y = 0;
	for (x = 0; x < width; x++) {
	    for (y = 0; y < height; y++) {
		if (mask[x][y])
		    signedDistance[x][y] = -Float.MAX_VALUE;
		else
		    signedDistance[x][y] = Float.MAX_VALUE;
	    }

	}
	this.signedDistance = ExactDistanceTransform(signedDistance);
	// renormalize(5);
    }

    /**
     * If we already have a signed-distance available, use it
     * 
     * @param signedDistance
     */
    public GenericLevelSet(float[][] signedDistance) {
	this.signedDistance = signedDistance;
	this.width = signedDistance.length;
	this.height = signedDistance[0].length;
    }

    public static float[][] maskToSignedDistance(boolean[][] mask) {
	int width = mask.length;
	int height = mask[0].length;
	float[][] grid = new float[width][height];
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		if (mask[x][y])
		    grid[x][y] = -Float.MAX_VALUE;
		else
		    grid[x][y] = Float.MAX_VALUE;
	    }

	}
	return (ExactDistanceTransform(grid));
    }

    public float[][] getSignedDistance() {
	return signedDistance;
    }

    public GenericLevelSet getComplement() {
	float[][] sdistance = new float[width][height];
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		sdistance[x][y] = -this.signedDistance[x][y];
	    }
	}
	return new GenericLevelSet(sdistance);
    }

    /**
     * Exact distance transformation in 1D taken from the Fiji project
     * 
     * @param f
     * @return
     */
    public static float[] EDTTransform1D(float[] f) {
	float d[] = new float[f.length];
	float[] fNeg = new float[f.length + 1];
	float[] zNeg = new float[f.length + 1];
	int[] yNeg = new int[f.length + 1];
	float[] fPos = new float[f.length + 1];
	float[] zPos = new float[f.length + 1];
	int[] yPos = new int[f.length + 1];

	fNeg[0] = -Float.MAX_VALUE;
	yNeg[0] = -1;
	zNeg[0] = Float.MAX_VALUE;
	int kNeg = 0, kPos = 0;
	fPos[0] = Float.MAX_VALUE;
	yPos[0] = -1;
	zPos[0] = Float.MAX_VALUE;

	float fval, fvalNeg, fvalPos, s;

	for (int q = 0; q < f.length; q++) {
	    fval = f[q];
	    fvalNeg = fval < 0 ? fval : 0f;
	    s = ((fvalNeg - q * q) - (fNeg[kNeg] - yNeg[kNeg] * yNeg[kNeg]))
		    / -2 / (q - yNeg[kNeg]);
	    for (;;) {
		// calculate the intersection
		s = ((fvalNeg - q * q) - (fNeg[kNeg] - yNeg[kNeg] * yNeg[kNeg]))
			/ -2 / (q - yNeg[kNeg]);
		if (s > zNeg[kNeg])
		    break;
		if (--kNeg < 0)
		    break;
	    }

	    kNeg++;
	    yNeg[kNeg] = q;
	    fNeg[kNeg] = fvalNeg;
	    zNeg[kNeg] = s;
	    fvalPos = fval > 0 ? fval : 0f;
	    for (;;) {
		// calculate the intersection
		s = ((fvalPos + q * q) - (fPos[kPos] + yPos[kPos] * yPos[kPos]))
			/ 2 / (q - yPos[kPos]);
		if (s > zPos[kPos])
		    break;
		if (--kPos < 0)
		    break;
	    }
	    kPos++;
	    yPos[kPos] = q;
	    fPos[kPos] = fvalPos;
	    zPos[kPos] = s;
	}
	zNeg[++kNeg] = Float.MAX_VALUE;
	zPos[++kPos] = Float.MAX_VALUE;

	int iNeg = 0, iPos = 0;
	for (int q = 0; q < f.length; q++) {
	    while (zNeg[iNeg + 1] < q)
		iNeg++;
	    while (zPos[iPos + 1] < q)
		iPos++;

	    d[q] = f[q] < 0 ? -(q - yNeg[iNeg]) * (q - yNeg[iNeg]) + fNeg[iNeg]
		    : (q - yPos[iPos]) * (q - yPos[iPos]) + fPos[iPos];
	    // d[q] = d[q]<0? 0.5f - (float)Math.sqrt(-d[q]): -0.5f + (float)
	    // Math.sqrt(d[q]);
	}

	return d;

    }

    /**
     * Advect the signed distance according to speed field along the normals of
     * the zero-level set, using fast marching
     * 
     * @param spd
     * @param time
     * @return
     */
    public GenericLevelSet advect(SpeedField spd, float time) {

	/**
	 * Basically compute the new signed distance function and create a new
	 * level object
	 */

	float[][] signedDistance = new float[width][height];

	return new GenericLevelSet(signedDistance);

    }

    public void advectUniformSpeed(double speed, double time) {
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		this.signedDistance[x][y] -= speed * time;
	    }
	}
    }

    /**
     * Perform the exact distance transformation given a signed distance
     * 
     * @param signedDistance
     * @return
     */
    public static float[][] ExactDistanceTransform(float[][] grid) {
	// Restore the signed distance

	float[][] newgrid = new float[grid.length][grid[0].length];

	float[] c = new float[grid.length];
	float[] r = new float[grid[0].length];
	for (int x = 0; x < grid.length; x++) {
	    for (int y = 0; y < grid[0].length; y++) {
		r[y] = grid[x][y] < 0 ? -Float.MAX_VALUE : Float.MAX_VALUE;
	    }
	    float[] d1 = EDTTransform1D(r);
	    for (int y = 0; y < grid[0].length; y++) {
		newgrid[x][y] = d1[y];
	    }
	}

	for (int y = 0; y < grid[0].length; y++) {
	    for (int x = 0; x < grid.length; x++) {
		c[x] = newgrid[x][y];
	    }
	    float[] d2 = EDTTransform1D(c);
	    for (int x = 0; x < grid.length; x++) {
		newgrid[x][y] = d2[x];
	    }
	}

	for (int x = 0; x < grid.length; x++) {
	    for (int y = 0; y < grid[0].length; y++) {
		newgrid[x][y] = (float) (newgrid[x][y] < 0 ? 0.5f - Math
			.sqrt(-newgrid[x][y]) : -0.5f
			+ Math.sqrt(newgrid[x][y]));
	    }
	}

	return newgrid;

    }

    public float get(int x, int y) {

	return this.signedDistance[Math.max(0, Math.min(x, width - 1))][Math
		.max(0, Math.min(y, height - 1))];
    }

    public float get(double x, double y) {
	// Interpolate a point

	int xfloor = (int) Math.floor(x);
	int xceil = xfloor+1;
	int yfloor = (int) Math.floor(y);
	int yceil = yfloor+1;

	return (float) ((1.0 * get(xfloor, yfloor) * (xceil - x) * (yceil - y)
		+ get(xceil, yfloor) * (x - xfloor) * (yceil - y)
		+ get(xfloor, yceil) * (xceil - x) * (y - yfloor) + get(xceil,
		yceil)
		* (x - xfloor) * (y - yfloor)))
		/ (xceil - xfloor) / (yceil - yfloor);
    }

    /**
     * Get central curvature
     * 
     * @param x
     * @param y
     * @return
     */
    public double getCurvature(int x, int y) {
	double dphiX = (get(x + 1, y) - get(x - 1, y)) / 2;
	double dphiY = (get(x, y + 1) - get(x, y - 1)) / 2;
	double dphiXX = ((get(x + 1, y) + get(x - 1, y) - (2 * get(x, y)))) / 4;
	double dphiYY = ((get(x, y + 1) + get(x, y - 1) - (2 * get(x, y)))) / 4;
	double dphiXY = (get(x + 1, y + 1) - get(x + 1, y - 1)
		- get(x - 1, y + 1) + get(x - 1, y - 1)) / 4;

	double curvature = (dphiXX * dphiY * dphiY + dphiYY * dphiX * dphiX - 2
		* dphiX * dphiY * dphiXY)
		/ (Math.pow(dphiX * dphiX + dphiY * dphiY, 1.5) + EPS);

	return curvature; // * deltaPhi;
    }

    /**
     * Central differences gradient at point (x,y)
     * 
     * @param x
     * @param y
     * @return
     */
    public float[] getGradient(int x, int y) {
	float dphiX = (get(x + 1, y) - get(x - 1, y)) / 2;
	float dphiY = (get(x, y + 1) - get(x, y - 1)) / 2;
	return new float[] { dphiX, dphiY };
    }

    public double[] getNormalVector(int x, int y) {
	return null;

    }

    public ArrayList<int[]> getZeroLevelSet() {
	ArrayList<int[]> points = new ArrayList<int[]>();
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		if (Math.abs(get(x, y)) <= 1) {
		    int[] p = { x, y };
		    points.add(p);
		}
	    }
	}
	return points;

    }

    public float l2Norm(float[] vector) {
	float sum = 0;
	for (int i = 0; i < vector.length; i++) {
	    sum += vector[i] * vector[i];
	}
	return (float) Math.sqrt(sum);

    }

    /**
     * Reinitialize the signed-distance function by solving the Eikonal equation
     * (slow)
     * 
     * @param reinitializations
     *            The number of iterations of reinitialization
     */
    public void renormalize(int reinitializations) {
	float maxVal = (float) Math.sqrt(width * width + height * height);

	float[][] newGrid;

	for (int i = 0; i < reinitializations; i++) {
	    newGrid = new float[width][height];

	    for (int x = 0; x < width; x++) {
		newGrid[x][0] = maxVal;
		newGrid[x][height - 1] = maxVal;
	    }
	    for (int y = 0; y < height; y++) {
		newGrid[0][y] = maxVal;
		newGrid[width - 1][y] = maxVal;
	    }

	    for (int x = 1; x < width - 1; x++)
		for (int y = 1; y < height - 1; y++) {
		    float xy = get(x, y);
		    float phiXPlus = (get(x + 1, y) - xy);
		    float phiXMinus = (xy - get(x - 1, y));
		    float phiYPlus = (get(x, y + 1) - xy);
		    float phiYMinus = (xy - get(x, y - 1));

		    float dXSquared = 0;
		    float dYSquared = 0;

		    // Please notice that the sign of "Sign" is the same as that
		    // of "signedDistance"
		    // Please consult the report for further information
		    if (signedDistance[x][y] > 0) {
			float max = Math.max(phiXMinus, 0);
			float min = Math.min(phiXPlus, 0);
			dXSquared = Math.max(max * max, min * min);

			max = Math.max(phiYMinus, 0);
			min = Math.min(phiYPlus, 0);
			dYSquared = Math.max(max * max, min * min);

		    } else {
			float max = Math.max(phiXPlus, 0);
			float min = Math.min(phiXMinus, 0);
			dXSquared = Math.max(max * max, min * min);

			max = Math.max(phiYPlus, 0);
			min = Math.min(phiYMinus, 0);
			dYSquared = Math.max(max * max, min * min);

		    }

		    float normSquared = dXSquared + dYSquared;

		    float norm = (float) Math.sqrt(normSquared);

		    float sign = (float) (signedDistance[x][y] / Math
			    .sqrt(signedDistance[x][y] * signedDistance[x][y]
				    + normSquared));

		    float t = (float) 0.3; // A stable CFL condition
		    newGrid[x][y] = signedDistance[x][y] - sign * (norm - 1)
			    * t;
		}
	    signedDistance = newGrid;
	}
    }

    /**
     * Set the value of a signedDistance point
     * 
     * @param x
     * @param y
     * @param val
     */
    public void set(int x, int y, float val) {
	this.signedDistance[x][y] = val;

    }

    /**
     * How much of the inner region in signedDistance 1 is inner in
     * signedDistance 2
     * 
     * @param grid1
     * @param grid2
     * @return Return the overlap
     */
    public static float overlap(float[][] grid1, float[][] grid2) {
	int w1 = grid1.length;
	int w2 = grid2.length;
	int h1 = grid1[0].length;
	int h2 = grid2[0].length;
	if (w1 != w2 || h1 != h2)
	    return 0;
	float area = 0;
	long over = 0;
	long in = 0;
	for (int x = 0; x < w1; x++) {
	    for (int y = 0; y < h1; y++) {
		if (grid1[x][y] <= 0) {
		    in++;
		    if (grid2[x][y] <= 0) {
			over++;
		    }

		}
	    }
	}
	area = (float) over / in;
	return area;
    }

    /**
     * Reinitialize the signed distance function to a mask
     * 
     * @param mask
     */
    public void reinitialize(boolean[][] mask) {
	this.width = mask.length;
	this.height = mask[0].length;
	this.signedDistance = new float[width][height];
	int x = 0, y = 0;
	for (x = 0; x < width; x++) {
	    for (y = 0; y < height; y++) {
		if (mask[x][y])
		    signedDistance[x][y] = -Float.MAX_VALUE;
		else
		    signedDistance[x][y] = Float.MAX_VALUE;
	    }

	}
	this.signedDistance = ExactDistanceTransform(signedDistance);
	// renormalize(5);
    }

    public float deltafunction(float x) {
	return (float) (Math.pow(Math.PI * (x * x + 1), -1));
    }

    public boolean[][] getMask() {
	boolean[][] mask = new boolean[width][height];
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		mask[x][y] = signedDistance[x][y] <= 0 ? true : false;
	    }
	}
	return mask;
    }

    public GenericLevelSet clone() {
	return this;
    }



}