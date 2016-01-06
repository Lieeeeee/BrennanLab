package ucla.brennanlab.imagej.util.levelsets;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ucla.brennanlab.imagej.util.stochastic.ImageStats;
import ucla.brennanlab.imagej.util.stochastic.SpeedField;

import java.awt.*;
import java.text.DecimalFormat;

public class RoiImageLevelSet extends GenericLevelSet {

    static final double MU = 4;
    public Overlay zeroset;
    public ImageStats imageStats;

    public float[][] previousgrid;
    ImageProcessor ip;
    Roi userDrawn;
    SpeedField sf;
    double priorMeanSpeed;
    double priorSpeedSD;
    private float[][] initialGuess; // The initial segmentation
    private float[][] probgrid; // Grid of prior probabilities

    /**
     * @param ip
     * @param priorprob
     */
    public RoiImageLevelSet(ImageProcessor ip, float[][] priorprob) {
        super(thresholdToMask(ip));
        // Construct without initial ROI
        this.ip = ip;
        width = ip.getWidth();
        height = ip.getHeight();

        this.setPriorInformation(0, Float.POSITIVE_INFINITY, null);
        updateImageStats();

    }

    public RoiImageLevelSet(Roi roi, ImageProcessor ip) {
        super(roiToMask(roi, ip.getWidth(), ip.getHeight()));
        this.ip = ip;
        width = ip.getWidth();
        height = ip.getHeight();
        this.userDrawn = roi;
        initialGuess = signedDistance;
        this.setPriorInformation(0, Float.POSITIVE_INFINITY, null);
        updateImageStats();

    }

    public RoiImageLevelSet(PolygonRoi roi, ImageProcessor ip) {
        super(roiToMask(roi, ip.getWidth(), ip.getHeight()));
        this.ip = ip;
        width = ip.getWidth();
        height = ip.getHeight();
        this.userDrawn = roi;
        initialGuess = signedDistance;
        this.setPriorInformation(0, Float.POSITIVE_INFINITY, null);
        updateImageStats();

    }

    public RoiImageLevelSet(boolean[][] mask) {
        super(mask);
    }

    /**
     * @param roi
     * @param width
     * @param height
     */
    public RoiImageLevelSet(Roi roi, int width, int height) {
        super(roiToMask(roi, width, height));

        initialGuess = signedDistance;

    }

    public static boolean[][] thresholdToMask(ImageProcessor ip) {
        boolean[][] mask = new boolean[ip.getWidth()][ip.getHeight()];
        ImageProcessor toBlur = ip.duplicate();
        GaussianBlur gb = new GaussianBlur();
        gb.blurGaussian(toBlur, 3,3,0.01);
        int x = 0, y = 0;
        for (x = 0; x < ip.getWidth(); x++) {
            for (y = 0; y < ip.getHeight(); y++) {
                if (toBlur.getPixelValue(x, y) < 5)
                    mask[x][y] = false;
                else
                    mask[x][y] = true;
            }

        }
        return mask;
    }

    public static boolean[][] roiToMask(Roi r, int width, int height) {
        boolean[][] mask = new boolean[width][height];
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (r.contains(x, y))
                    mask[x][y] = true;
                else
                    mask[x][y] = false;
            }
        }
        return mask;

    }

    static public RoiImageLevelSet genericToRoi(GenericLevelSet ls) {
        boolean[][] mask = new boolean[ls.width][ls.height];
        for (int x = 0; x < ls.width; x++) {
            for (int y = 0; y < ls.height; y++) {
                mask[x][y] = ls.get(x, y) <= 0 ? true : false;
            }
        }
        return new RoiImageLevelSet(mask);
    }

    /**
     * Save evenly spaced X-Y coordinates of the current level set
     *
     * @return
     * @url{http://rsb.info.nih.gov/ij/plugins/download/Path_Writer.java
     */
    public static int[][] returnPath(Roi roi) {

        int n = ((PolygonRoi) roi).getNCoordinates();
        int[] x = ((PolygonRoi) roi).getXCoordinates();
        int[] y = ((PolygonRoi) roi).getYCoordinates();
        Rectangle r = roi.getBounds();
        int xbase = r.x;
        int ybase = r.y;
        boolean areaPath = roi.getType() <= Roi.TRACED_ROI;
        double length = 0.0;
        double segmentLength;
        int xdelta, ydelta;
        double[] segmentLengths = new double[n];
        int[] dx = new int[n];
        int[] dy = new int[n];
        for (int i = 0; i < (n - 1); i++) {
            xdelta = x[i + 1] - x[i];
            ydelta = y[i + 1] - y[i];
            segmentLength = Math.sqrt(xdelta * xdelta + ydelta * ydelta);
            length += segmentLength;
            segmentLengths[i] = segmentLength;
            dx[i] = xdelta;
            dy[i] = ydelta;
        }
        if (areaPath) {
            xdelta = x[0] - x[n - 1];
            ydelta = y[0] - y[n - 1];
            segmentLength = Math.sqrt(xdelta * xdelta + ydelta * ydelta);
            length += segmentLength;
            segmentLengths[n - 1] = segmentLength;
            dx[n - 1] = xdelta;
            dy[n - 1] = ydelta;
        }

        int size = (int) (1.1 * length);
        int pathLength = 0;
        double[] xpath = new double[size];
        double[] ypath = new double[size];
        double leftOver = 1.0;
        double distance = 0.0;
        int index = -1;
        for (int i = 0; i < n; i++) {
            double len = segmentLengths[i];
            if (len == 0.0)
                continue;
            double xinc = dx[i] / len;
            double yinc = dy[i] / len;
            double start = 1.0 - leftOver;
            double rx = xbase + x[i] + start * xinc;
            double ry = ybase + y[i] + start * yinc;
            double len2 = len - start;
            int n2 = (int) len2;

            for (int j = 0; j <= n2; j++) {
                index++;
                if (index < xpath.length) {
                    xpath[index] = rx;
                    ypath[index] = ry;
                    pathLength = index + 1;
                }
                rx += xinc;
                ry += yinc;
            }
            distance += len;
            leftOver = len2 - n2;
        }

        /**
         * OK, but we want to also exclude boundary points from our path, so
         * let's detect where we find boundary points, excise them, and redo the
         * ordering of the remaining points as needed
         *
         */

        int[][] path = new int[pathLength][2];

        for (int i = 0; i < pathLength; i++) {

            int x1 = (int) Math.round(xpath[i]);
            int y1 = (int) Math.round(ypath[i]);
            path[i][0] = x1 - 1;
            path[i][1] = y1 - 1;
        }

        return (path);

    }

    /**
     * Calculates the energy and updates the stats for a given segmentation
     *
     * @param grid
     * @param MU
     * @return
     */
    public ImageStats calculateEnergy(float[][] grid, double MU) {
        ImageStats s = new ImageStats();
        /*************************************************
         * Find conditional pixel intensity distributions
         ************************************************/

        double im = 0;
        double iim = 0;
        double om = 0;
        double iv = 0;
        double iiv = 0;
        double ov = 0;

        int N = 0, M = 0, x = 0, y = 0, L = 0;
        for (x = 0; x < width; x++) {
            for (y = 0; y < height; y++) {
                double pixval = ip.getPixelValue(x, y);
                double gridval = this.get(x, y);
                if (gridval > 0) {
                    om += pixval;
                    ov += pixval * pixval;
                    N++;
                } else {
                    M++;
                    im += pixval;
                    iv += pixval * pixval;
                    if (this.get(x, y) > -20) {
                        iim += pixval;
                        iiv += pixval * pixval;
                        L++;
                    }

                }
            }
        }

        s.innerarea = 1.0 * M / width / height;
        s.innermean = im / M;
        s.outermean = om / N;
        s.innerinnermean = iim / L;
        s.innersd = Math.sqrt(iv / M - Math.pow(im / M, 2));
        s.innerinnersd = Math.sqrt(iiv / L - Math.pow(iim / L, 2));
        s.outersd = Math.sqrt(ov / N - Math.pow(om / N, 2));

        s.energy = 0;
        double length = 0;
        for (x = 0; x < width; x++) {
            for (y = 0; y < height; y++) {
                s.energy += 0.5
                        * (grid[x][y] < 0 ? Math.log(2 * Math.PI
                        * s.innerinnersd * s.innerinnersd + EPS)
                        + Math
                        .pow(
                                (ip.getPixelValue(x, y) - s.innerinnermean)
                                        / s.innerinnersd, 2)
                        - 2 * Math.log(probgrid[x][y] + EPS)
                        : Math.log(2 * Math.PI * s.outersd * s.outersd
                        + EPS)
                        + Math
                        .pow(
                                (ip.getPixelValue(x, y) - s.outermean)
                                        / s.outersd, 2))
                        - 2 * Math.log(1 - probgrid[x][y] + EPS);
                length += Math.abs(grid[x][y]) < 1 ? l2Norm(getGradient(x, y))
                        : 0;
            }
        }
        s.energy += MU * length;
        return s;

    }

    public double calculateOdds(float[][] segmentation1,
                                float[][] segmentation2, ImageProcessor ip) {
        ImageStats s1 = calculateEnergy(segmentation1, MU);
        ImageStats s2 = calculateEnergy(segmentation2, MU);
        return Math.pow(1 + Math.exp(s1.energy - s2.energy), -1);
    }

    public ImagePlus getCurvatureImage() {
        ImageProcessor ip2 = new FloatProcessor(width, height);
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                ip2.putPixelValue(x, y, getCurvature(x, y));
            }
        }
        ip2.resetMinAndMax();
        ImagePlus imp = new ImagePlus("Curvature", ip2);
        imp.resetDisplayRange();
        return imp;
    }

    public FloatProcessor getLSFloatProcessor() {
        FloatProcessor ip2 = new FloatProcessor(width, height);
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                ip2.setf(x, y, signedDistance[x][y]);
            }
        }
        ip2.resetMinAndMax();
        ip2.crop();
        return (ip2);
    }

    public ImagePlus getLSImagePlus() {
        FloatProcessor ip2 = getLSFloatProcessor();
        ImagePlus imp = new ImagePlus("Signed distance", ip2);
        return imp;

    }

    /**
     * Get the ROI of the current level set
     *
     * @param reset Reset the signed distance function after obtaining the ROI
     * @return
     */
    public PolygonRoi getRoi(boolean reset) {
        // Returns the ROI corresponding to the zero level set
        ImageProcessor ip2 = this.getLSFloatProcessor();
        int minpt[] = {0, 0};
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (signedDistance[x][y] < signedDistance[minpt[0]][minpt[1]]) {
                    minpt[0] = x;
                    minpt[1] = y;
                }
            }
        }

        Wand w = new Wand(ip2);
        w.autoOutline(minpt[0], minpt[1], signedDistance[minpt[0]][minpt[1]],
                0, Wand.EIGHT_CONNECTED);
        if (w.npoints == 0) {
            return null;
        }

        Polygon p = new Polygon(w.xpoints, w.ypoints, w.npoints);
        PolygonRoi poly = new PolygonRoi(p, Roi.POLYGON);
        if (reset) {
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    if (poly.contains(x, y))
                        signedDistance[x][y] = -Float.MAX_VALUE;
                    else
                        signedDistance[x][y] = Float.MAX_VALUE;
                }

            }

            this.signedDistance = ExactDistanceTransform(this.signedDistance);
            renormalize(2);
        }

        return poly;
    }

    public void recalculateLevelSets(PolygonRoi roi) {
        int x = 0, y = 0;
        for (x = 0; x < width; x++) {
            for (y = 0; y < height; y++) {
                if (roi.contains(x, y))
                    signedDistance[x][y] = -Float.MAX_VALUE;
                else
                    signedDistance[x][y] = Float.MAX_VALUE;
            }

        }
        this.previousgrid = signedDistance;
        this.signedDistance = ExactDistanceTransform(this.signedDistance);
        renormalize(2);
        updateImageStats();
    }

    public void reloadImage(ImageProcessor ip2) {
        this.ip = ip2;
        this.previousgrid = this.signedDistance;
    }

    /**
     * Set the prior information for energy minimization
     *
     * @param priorspeed
     * @param priorspeedsd
     * @param probgrid
     */
    public void setPriorInformation(double priorspeed, double priorspeedsd,
                                    float[][] probgrid) {
        // do pre-conversion to pixels before calling this function
        this.priorMeanSpeed = priorspeed;
        this.priorSpeedSD = priorspeedsd;
        if (probgrid != null)
            this.probgrid = probgrid;
        else {
            this.probgrid = new float[width][height];
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    this.probgrid[x][y] = (float) 0.5;
                }
            }
        }
    }

    public float[][] stepMAP(int steps, float MU) {
        return null;
    }

    /**
     * Minimize the image energy
     *
     * @param steps   number of steps to run gradient descent
     * @param verbose Whether or not to output a lot
     */
    public void solveMAP(int steps, boolean verbose, double MU) {
        // Solves for the Maximum A-Posterior estimate for the level set using
        // gradient descent
        double timegrain = .7;
        double[] dT = new double[7];
        int time = 0;
        double[] energies = new double[steps];
        while (time++ < steps) {
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    double dTemp = timegrain
                            / (1 + signedDistance[x][y] * signedDistance[x][y]);
                    double dPhiX = get(x + 1, y) - get(x, y);
                    double dPhiY = 0.5 * (get(x, y - 1) - get(x, y + 1));
                    dT[0] = 1.0 / Math
                            .sqrt(EPS + dPhiX * dPhiX + dPhiY * dPhiY);
                    dPhiX = get(x, y) - get(x - 1, y);
                    dPhiY = (get(x - 1, y + 1) - get(x - 1, y - 1)) * 0.5;
                    dT[1] = 1.0 / Math
                            .sqrt(EPS + dPhiX * dPhiX + dPhiY * dPhiY);
                    dPhiX = (get(x + 1, y) - get(x - 1, y)) * 0.5;
                    dPhiY = get(x, y + 1) - get(x, y);
                    dT[2] = 1 / Math.sqrt(EPS + dPhiX * dPhiX + dPhiY * dPhiY);
                    dPhiX = (get(x + 1, y - 1) - get(x - 1, y - 1)) * 0.5;
                    dPhiY = get(x, y) - get(x, y - 1);//
                    dT[3] = 1 / Math.sqrt(EPS + dPhiX * dPhiX + dPhiY * dPhiY);
                    dT[4] = 0.5
                            * Math
                            .pow(
                                    (ip.getPixelValue(x, y) - imageStats.outermean)
                                            / (EPS + imageStats.outersd
                                            * imageStats.outersd),
                                    2)
                            + 0.5
                            * Math.log(2 * Math.PI * imageStats.outersd
                            * imageStats.outersd + EPS);
                    dT[5] = 0.5
                            * Math
                            .pow(
                                    (ip.getPixelValue(x, y) - imageStats.innerinnermean)
                                            / (EPS + imageStats.innerinnersd
                                            * imageStats.innerinnersd),
                                    2)
                            + 0.5
                            * Math.log(2 * Math.PI * imageStats.innerinnersd
                            * imageStats.innerinnersd + EPS);
                    dT[6] = -Math.log(probgrid[x][y] + EPS)
                            + Math.log(1 - probgrid[x][y] + EPS);
                    signedDistance[x][y] = (float) ((get(x, y) + dTemp
                            * (MU
                            * (get(x + 1, y) * dT[0] + get(x - 1, y)
                            * dT[1] + get(x, y + 1) * dT[2] + get(
                            x, y - 1)
                            * dT[3]) - dT[4] + dT[5] + dT[6])) / (1.0 + dTemp
                            * MU * (dT[0] + dT[1] + dT[2] + dT[3])));
                    ;
                }

            }
            for (int i = 1; i < width - 1; i++) {
                signedDistance[i][0] = signedDistance[i][1];
                signedDistance[i][height - 1] = signedDistance[i][height - 2];
            }
            for (int j = 1; j < height - 1; j++) {
                signedDistance[0][j] = signedDistance[1][j];
                signedDistance[width - 1][j] = signedDistance[width - 2][j];
            }
            renormalize(1);
            updateImageStats();
            energies[time - 1] = imageStats.energy;

            IJ.showProgress((double) time / steps);

            if (imageStats.innerarea == 0 || imageStats.innerinnersd == 0
                    || imageStats.outersd == 0
                    || Double.isNaN(imageStats.energy)
                    || Double.isNaN(imageStats.innermean)) {// detect
                // collapse
                // IJ.showMessage("collapse " + innerarea);
                break;
            }
            DecimalFormat df = new DecimalFormat("#.###");
            if (verbose && time % 20 == 0) {
                // draw location of current ROI too?

                IJ.log("Slice " + ip.getSliceNumber() + " step "
                        + df.format(time) + "     area: "
                        + df.format(imageStats.innerarea * 100)
                        + "%\n    Energy: " + df.format(imageStats.energy)
                        + "\n    mean in "
                        + df.format(imageStats.innerinnermean) + " mean out "
                        + df.format(imageStats.outermean) + "\n    sd in "
                        + df.format(imageStats.innerinnersd) + " sd out "
                        + df.format(imageStats.outersd));
            }

        }
    }

    /**
     * Sample about MAP contour with gradient descent steps
     *
     * @param MU
     * @param steps
     * @return A level set a certain number of gradient descent steps away
     */
    public float[][] sampleMAP(float MU, float steps) {
        return initialGuess;

    }

    /**
     * Updates the class stat variables
     */
    public void updateImageStats() {
        ImageStats s = calculateEnergy(signedDistance, MU);
        this.imageStats = s;
    }

}
