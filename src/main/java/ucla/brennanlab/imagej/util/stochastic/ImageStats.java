package ucla.brennanlab.imagej.util.stochastic;

import ij.process.ImageProcessor;
import ucla.brennanlab.imagej.util.levelsets.GenericLevelSet;

/**
 * .
 * Compute the regional image stats and energy associated with a segmentation
 *
 * @author josh
 */
public class ImageStats {
    public double innermean, innerinnermean;
    public double outermean;
    public double innersd, innerinnersd;
    public double outersd;
    public double innerarea;
    public double outerarea;
    public double energy;
    public double priorinnermean = 0;
    public double pooledmean = 0;
    public double pooledsd = 0;
    public double priorinnersd = Double.MAX_VALUE;
    public double EPS = 1e-10;
    private double prioroutermean = 0;
    private double prioroutersd;


    /**
     * Initialize imagestats with prior information
     *
     * @param priormean
     * @param priorsd
     */

    /**
     * Non-informative prior
     */
    public ImageStats() {
        this(0, Double.MAX_VALUE, 0, Double.MAX_VALUE);
    }

    /**
     * Informative prior
     *
     * @param priormean
     * @param priorsd
     */
    public ImageStats(double priorinnermean, double priorinnersd,
                      double prioroutermean, double prioroutersd) {

        this.priorinnermean = priorinnermean;
        this.innermean = priorinnermean;
        this.innerinnermean = priorinnermean;
        this.innerinnersd = priorinnersd;
        this.innersd = priorinnersd;
        this.priorinnersd = priorinnersd;
        this.outermean = prioroutermean;
        this.prioroutermean = prioroutermean;
        this.prioroutersd = prioroutersd;
        this.outersd = prioroutersd;
        this.pooledmean = (priorinnermean + prioroutermean) / 2;
        this.pooledsd = Math.sqrt(priorinnersd * priorinnersd + prioroutersd * prioroutersd);

    }


    /**
     * Compute the energy of a particular mask
     *
     * @param pixelvals
     * @param mask
     * @param priors
     */
    public void calculateEnergy(float[][] pixelvals, boolean[][] mask,
                                float[][] priors) {

    }


    public void updateNarrowBandStats(ImageProcessor ip, boolean[][] mask, int band) {
        int width = ip.getWidth();
        float[][] signedDistance = GenericLevelSet.maskToSignedDistance(mask);
        int height = ip.getHeight();
        if (width != mask.length || height != mask[0].length)
            return;

        float sinner = 0, sinnerinner = 0, ssinner = 0, ssinnerinner = 0, souter = 0, ssouter = 0;
        int ninner = 0, ninnerinner = 0, nouter = 0;

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (mask[x][y]) {
                    if (Math.abs(signedDistance[x][y]) < band) {
                        ninnerinner++;
                        sinnerinner += ip.getPixelValue(x, y);
                        ssinnerinner += ip.getPixelValue(x, y) * ip.getPixelValue(x, y);
                    }
                    ninner++;
                    sinner += ip.getPixelValue(x, y);
                    ssinner += ip.getPixelValue(x, y) * ip.getPixelValue(x, y);
                } else {
                    nouter++;
                    souter += ip.getPixelValue(x, y);
                    ssouter += ip.getPixelValue(x, y) * ip.getPixelValue(x, y);
                }
            }
        }


        this.innerarea = (float) ninner / width / height;

        this.outerarea = 1 - this.innerarea;

        double imean = sinner / ninner;
        double iimean = sinnerinner / ninnerinner;
        double imse = ssinner / ninner - imean * imean;
        double iimse = ssinnerinner / ninnerinner - iimean * iimean;

        double omean = souter / nouter;
        double omse = ssouter / nouter - omean * omean;

        if (ninner == 0) {
            this.innersd = Double.MAX_VALUE;
            this.innermean = 0;
            this.outermean = (this.prioroutermean / this.prioroutersd + souter
                    / omse)
                    / (1 / this.prioroutersd / this.prioroutersd + nouter / omse);

            this.outersd = Math.pow(1 / this.prioroutersd / this.prioroutersd
                    + nouter / omse, -0.5);
            this.outersd = Math.pow(omse, -0.5);
            this.pooledmean = this.outermean;
            this.pooledsd = this.outersd;
            this.innerinnersd = Double.MAX_VALUE;
            this.innerinnermean = 0;
            return;
        } else if (nouter == 0) {
            this.outermean = 0;
            this.outersd = Double.MAX_VALUE;
            this.innermean = (this.priorinnermean / this.priorinnersd + sinner
                    / imse)
                    / (1 / (this.priorinnersd * this.priorinnersd + EPS) + ninner / imse);

            this.innersd = Math.pow(1 / this.priorinnersd / this.priorinnersd
                    + ninner / imse, -0.5);
            this.innersd = Math.pow(imse, -0.5);
            this.pooledmean = this.innermean;
            this.pooledsd = this.innersd;
            return;
        }


        this.innermean = (this.priorinnermean / this.priorinnersd + sinner
                / imse)
                / (1 / (this.priorinnersd * this.priorinnersd + EPS) + ninner / imse);
        this.innerinnermean = (this.priorinnermean / this.priorinnersd + sinnerinner
                / iimse)
                / (1 / (this.priorinnersd * this.priorinnersd + EPS) + ninnerinner / iimse);

        this.innersd = Math.pow(imse, 0.5);
        this.innerinnersd = Math.pow(iimse, 0.5);

        this.outermean = (this.prioroutermean / this.prioroutersd + souter
                / omse)
                / (1 / this.prioroutersd / this.prioroutersd + nouter / omse);

        this.outersd = Math.pow(1 / this.prioroutersd / this.prioroutersd
                + nouter / omse, 0.5);
        this.outersd = Math.pow(omse, 0.5);
        this.innersd = this.outersd;
        this.pooledmean = (nouter * this.outermean + ninner * this.innermean) / (nouter + ninner);
        this.pooledsd = (Math.sqrt(ninner * this.innersd * this.innersd + nouter * this.outersd * this.outersd)) / (ninner + nouter);
        return;


    }

    /**
     * Update the Gaussian Image stats according to a mask, weighted against
     * the prior information
     *
     * @param ip
     * @param mask
     */
    public void updateStats(ImageProcessor ip, boolean[][] mask) {
        updateNarrowBandStats(ip, mask, Integer.MAX_VALUE);
    }

    /**
     * Update stats according to signed distance function
     *
     * @param ip
     * @param signedDistance
     */
    public final void updateStats(ImageProcessor ip, float[][] signedDistance) {
        int width = signedDistance.length;
        int height = signedDistance[0].length;
        boolean[][] mask = new boolean[width][height];
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                mask[x][y] = signedDistance[x][y] < 0 ? true : false;
            }
        }
        updateStats(ip, mask);
    }

}