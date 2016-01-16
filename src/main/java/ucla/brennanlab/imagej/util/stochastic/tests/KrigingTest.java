package ucla.brennanlab.imagej.util.stochastic.tests;

/**
 * Use values in ROI to interpolate rest of image
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import org.ujmp.core.DenseMatrix;
import org.ujmp.core.Matrix;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;
import ucla.brennanlab.imagej.util.stochastic.Kriging2DLattice;

public class KrigingTest implements PlugInFilter {

    ImagePlus imp;
    String arg;
    Roi roi;

    public void run(ImageProcessor ip) {

        int gridwidth = 5;// reduce the dimensions down to 2x2 blocks
        int gridheight = 5;
        int bandwidth = 30; // band for interpolation, all outside is set to
        // mean

        IJ.log("Starting kriging test");

        int width = imp.getWidth();
        int height = imp.getHeight();

        int newWidth = (int) Math.ceil(1.0 * width / gridwidth);
        int newHeight = (int) Math.ceil(1.0 * height / gridheight);

        boolean[][] mask = new boolean[width][height];

        int location = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (roi.contains(x, y)) {
                    mask[x][y] = true;
                    location++;
                } else {
                    mask[x][y] = false;
                }
            }
        }

        ImplicitShape2D signedDistance = new ImplicitShape2D(mask);

        // Masks over the new coarsened grid

        boolean[][] predictMask = new boolean[newWidth][newHeight];
        boolean[][] interpMask = new boolean[newWidth][newHeight];
        int inside = 0, outside = 0;

        for (int y = 0; y < newHeight; y++) {
            for (int x = 0; x < newWidth; x++) {
                // location on original grid for point in center of new grid
                double originalx = (gridwidth + 1.0) / 2 + gridwidth * x;
                double originaly = (gridheight + 1.0) / 2 + gridheight * y;
                if (signedDistance.get(originalx, originaly) <= 0
                        && signedDistance.get(originalx, originaly) > -bandwidth) {
                    interpMask[x][y] = false;
                    predictMask[x][y] = true;
                    inside++;
                } else if (signedDistance.get(originalx, originaly) > 0
                        && signedDistance.get(originalx, originaly) < bandwidth) {
                    interpMask[x][y] = true;
                    predictMask[x][y] = false;
                    outside++;
                } else {
                    predictMask[x][y] = false;
                    interpMask[x][y] = false;
                }
            }
        }

        IJ.log("Predicting using " + inside + " points");
        IJ.log("Interpolating " + outside + " points");
        Matrix locations = DenseMatrix.Factory.zeros(inside, 2);
        Matrix response = DenseMatrix.Factory.zeros(inside, 1);
        Matrix covariates = DenseMatrix.Factory.zeros(inside, 1);
        Matrix interpolationLocations = DenseMatrix.Factory.zeros(outside, 2);
        Matrix interpolationCovariates = DenseMatrix.Factory.zeros(outside, 1);

        int locationpred = 0;
        int locationinterp = 0;
        for (int y = 0; y < newHeight; y++) {
            for (int x = 0; x < newWidth; x++) {
                if (predictMask[x][y]) {
                    locations.setAsInt(x, locationpred, 0);
                    locations.setAsInt(y, locationpred, 1);
                    covariates.setAsInt(1, locationpred, 0);
                    // compute response
                    double totalresponse = 0;

                    for (int yminor = 0; yminor < gridheight; yminor++) {
                        for (int xminor = 0; xminor < gridwidth; xminor++) {
                            totalresponse += ip.getPixelValue(
                                    Math.min(Math
                                                    .max(xminor + x * gridwidth, 0),
                                            width - 1), Math.min(Math.max(
                                            yminor + y * gridheight, 0),
                                            height - 1));
                        }
                    }
                    response.setAsDouble(
                            totalresponse / gridheight / gridwidth,
                            locationpred, 0);
                    locationpred++;

                } else if (interpMask[x][y]) {
                    interpolationLocations.setAsInt(x, locationinterp, 0);
                    interpolationLocations.setAsInt(y, locationinterp, 1);
                    interpolationCovariates.setAsInt(1, locationinterp, 0);
                    locationinterp++;
                }
            }
        }


        Matrix betaprior = DenseMatrix.Factory.zeros(1, 1);
        Matrix V = DenseMatrix.Factory.ones(1, 1);
        V.setAsDouble(Double.MIN_VALUE, 0, 0);

        long startTime = System.nanoTime();
        long endTime;
        Kriging2DLattice krig;
        krig = new Kriging2DLattice(betaprior, V,3,3);

        startTime = System.nanoTime();
        try {

            // perform kriging and output the results
            IJ.log("Starting kriging to estimate beta");
            IJ.showStatus("Kriging to estimate beta...");
            krig.computeBeta();
            IJ.showStatus("Done kriging to estimate");

        } finally {
            endTime = System.nanoTime();

        }

        IJ.log("Finished kriging to estimate beta =" + krig.getBeta() + "+-"
                + Math.sqrt(krig.getBetaVariance().getAsDouble(0, 0)) + " in "
                + (float) (endTime - startTime) / 1000000000 + " s");

        startTime = System.nanoTime();

        try {

            IJ.showStatus("interpolation...");
            krig.interpolate(interpolationLocations, interpolationCovariates);
        } finally {
            endTime = System.nanoTime();
            IJ.log("Finished interpolating in " + (float) (endTime - startTime)
                    / 1000000000 + " s");
        }
        IJ.showStatus("Done interpolating");

        Matrix interpolatedMean = krig.getPredicted();
        ImageProcessor ip2 = ip.createProcessor(width, height);

        // Set values for new processor

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (signedDistance.get(x, y) < 1) {
                    ip2.setf(x, y, ip.getPixelValue(x, y));
                } else {
                    ip2.setf(x, y, krig.getBeta().floatValue());
                }

            }
        }

        // plug in the interpolated values

        Matrix interpolationMatrix = DenseMatrix.Factory.zeros(newWidth,
                newHeight);
        // interpolationMatrix.showGUI();
        double mean = krig.getBeta().doubleValue();
        for (int y = 0; y < newHeight; y++) {
            for (int x = 0; x < newWidth; x++) {
                interpolationMatrix.setAsDouble(mean, x, y);
            }
        }

        for (int j = 0; j < interpolationLocations.getRowCount(); j++) {
            int gridlocx = interpolationLocations.getAsInt(j, 0);
            int gridlocy = interpolationLocations.getAsInt(j, 1);
            interpolationMatrix.setAsDouble(interpolatedMean.getAsDouble(j, 0),
                    gridlocx, gridlocy);
        }

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (signedDistance.get(x, y) > 0) {
                    // interpolate using the interpolationMatrix
                    int gridlocx = (int) Math.floor(1.0 * x / gridwidth);
                    int gridlocy = (int) Math.floor(1.0 * y / gridheight);
                    if (x > width - gridwidth / 2 - 3 || x < gridwidth / 2 + 1 || y < gridheight / 2 + 1 || y > height - gridheight / 2 - 3) {
                        // just set equal to the value within the coarse grid
                        ip2.setf(x, y, interpolationMatrix.getAsFloat(gridlocx,
                                gridlocy));

                    } else {
                        // want to interpolate using the coarse grid
                        // need to determine which quadrant we are in
                        // to set up the corners to use for interpolation

                        int xoffset = x % gridwidth;
                        int yoffset = y % gridheight;

                        int downx, upx, downy, upy;

                        double f[] = new double[4];
                        if (xoffset <= gridwidth / 2) {
                            // on the left
                            if (yoffset <= gridheight / 2) {
                                // bottom left
                                upx = gridlocx;
                                upy = gridlocy;
                                downx = Math.max(0, gridlocx - 1);
                                downy = Math.max(0, gridlocy - 1);

                            } else {
                                upx = Math.min(gridlocx + 1, newWidth - 1);
                                downx = gridlocx;
                                upy = gridlocy;
                                downy = Math.max(0, gridlocy - 1);
                                // top left

                            }
                        } else {
                            // on the right
                            if (yoffset <= gridheight / 2) {
                                // bottom right

                                downx = gridlocx;
                                upx = Math.min(newWidth - 1, gridlocx + 1);
                                upy = gridlocy;
                                downy = Math.max(0, gridlocy - 1);
                            } else {

                                // top right
                                // (x,x+1) (y,y+1)
                                downx = gridlocx;
                                upx = Math.min(gridlocx + 1, newWidth - 1);
                                downy = gridlocy;
                                upy = Math.min(gridlocy + 1, newHeight - 1);

                            }

                        }
                        f[0] = interpolationMatrix.getAsFloat(downx, downy);
                        f[1] = interpolationMatrix.getAsFloat(downx, upy);
                        f[2] = interpolationMatrix.getAsFloat(upx, downy);
                        f[3] = interpolationMatrix.getAsFloat(upx, upy);
                        float value = (float) bilinearinterp(new double[]{
                                        downx, upx}, new double[]{downy, upy}, f,
                                new double[]{(x - 1.0) / gridwidth, (y - 1.0) / gridheight});
                        ip2.setf(x, y, Float.isNaN(value) || Float.isInfinite(value) ? krig.getBeta().floatValue() : value);

                    }
                }
            }
        }

        ImagePlus interpolatedImage = new ImagePlus("Interpolated", ip2);
        interpolatedImage.show();
        interpolatedImage.updateAndDraw();

        IJ.log("Done with kriging test");

        // krig.getPredictionRootVariance().showGUI();
        krig = null;

        return;
    }

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.arg = arg;
        this.roi = imp.getRoi();
        return DOES_32 + ROI_REQUIRED;
    }

    // [x1 , x2] [y1, y2] [f11, f12, f21, f22] [x y]
    private double bilinearinterp(double[] x, double[] y, double[] f,
                                  double[] pt) {
        // return f[0];
        double ans = f[0] * (x[1] - pt[0]) * (y[1] - pt[1]) + f[2]
                * (pt[0] - x[0]) * (y[1] - pt[1]) + f[1] * (x[1] - pt[0])
                * (pt[1] - y[0]) + f[3] * (pt[0] - x[0]) * (pt[1] - y[0]);
        ans /= (x[1] - x[0]) * (y[1] - y[0]);
        return ans;

    }

}
