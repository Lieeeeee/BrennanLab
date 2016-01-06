package ucla.brennanlab.imagej.util.stochastic;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.Opener;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import javax.swing.*;

public class Gibbs_MRF_sampler implements PlugInFilter {

    static int iters = 100;
    protected ImagePlus imp;
    protected String arg;
    double theta = 1;
    double psi = 1;
    double phi = 1;
    int width;
    int height;
    String trainingfile;

    private double Energy(double val, double[][] Y, int i, int j, double[][] X) { // Gives
        // the
        // energy
        // associated
        // with changing the value of a single node i,j to a new value value
        double energy;
        if (i == 0) {
            if (j == 0) {

                energy = (phi(val, Y[i + 1][j]) + phi(val, Y[i][j + 1]));
            } else if (j == (Y[0].length - 1))
                energy = (phi(val, Y[i + 1][j]) + phi(val, Y[i][j - 1]));
            else
                energy = (phi(val, Y[i + 1][j]) + phi(val, Y[i][j - 1]) + phi(
                        val, Y[i][j + 1]));
        } else if (i == Y.length - 1) {
            if (j == 0) {
                energy = (phi(val, Y[i - 1][j]) + phi(val, Y[i][j + 1]));
            } else if (j == Y[0].length - 1) {
                energy = (phi(val, Y[i - 1][j]) + phi(val, Y[i][j - 1]));
            } else {
                energy = (phi(val, Y[i - 1][j]) + phi(val, Y[i][j - 1]) + phi(
                        val, Y[i][j + 1]));
            }
        } else {
            if (j == 0) {
                energy = (phi(val, Y[i - 1][j]) + phi(val, Y[i + 1][j]) + phi(
                        val, Y[i][j + 1]));
            } else if (j == Y[0].length - 1) {
                energy = (phi(val, Y[i - 1][j]) + phi(val, Y[i + 1][j]) + phi(
                        val, Y[i][j - 1]));
            } else {
                energy = (phi(val, Y[i - 1][j]) + phi(val, Y[i + 1][j]) + phi(
                        val, Y[i][j - 1] + Y[i][j + 1]));
            }
        }
        energy += psi(val, X[i][j]) + theta(val);
        return (energy);

    }

    private double[][] getMatrix(ImageProcessor ip) {
        int width = ip.getWidth();
        int height = ip.getHeight();
        double[][] matrix = new double[width][height];
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                matrix[x][y] = ip.getPixel(x, y);
            }
        }

        return matrix;
    }

    private double[] isingprob(double[][] Y, int i, int j, double[][] X) {
        double up = Math.exp(-Energy(1, Y, i, j, X));
        double down = Math.exp(-Energy(-1, Y, i, j, X));
        double zero = Math.exp(-Energy(0, Y, i, j, X));
        double Z = up + down + zero;
        double prob[] = {up / Z, zero / Z, down / Z};
        return prob;
    }

    private void learn(double[][] Y, double[][] X) {

        double thetahat = 0;
        double psihat = 0;
        double phihat = 0;
        int edgecount = 0;
        for (int i = 0; i < Y.length - 1; i++) { // sum over left interior nodes
            // and edges
            for (int j = 0; j < Y[0].length - 1; j++) {

                thetahat += theta(Y[i][j]);
                psihat += psi(Y[i][j], X[i][j]);
                phihat += phi(Y[i][j], Y[i + 1][j])
                        + phi(Y[i][j], Y[i + 1][j + 1])
                        + phi(Y[i][j], Y[i][j + 1]);
                edgecount += 3;
            }
            thetahat += theta(Y[i][Y[0].length - 1]);
            phihat += phi(Y[i][Y[0].length - 1], Y[i + 1][Y[0].length - 1]); // right
            // edge
            edgecount += 2;
        }
        for (int j = 0; j < Y[0].length - 1; j++) {
            phihat += phi(Y[Y.length - 1][j], Y[Y.length - 1][j + 1]);
            edgecount += 2;
        }
        this.psi = psihat / Y.length / Y[0].length;
        this.phi = phihat / edgecount;
    }

    public FloatProcessor MatrixToFloatProcessor(double[][] matrix) {
        float[][] pixels = new float[matrix.length][matrix[1].length];

        for (int x = 0; x < matrix.length; x++)
            for (int y = 0; y < matrix[1].length; y++)
                pixels[x][y] = (float) matrix[x][y];

        FloatProcessor fp = new FloatProcessor(pixels);
        return fp;
    }

    protected GenericDialog OptionSelector() {
        GenericDialog gd = new GenericDialog("Options", IJ.getInstance());
        gd.addNumericField("Gibbs iterations:", iters, 0);
        gd.addNumericField("phi parameter:", phi, 2);
        gd.addNumericField("theta parameter:", theta, 2);
        gd.addNumericField("psi parameter:", psi, 2);
        gd.addCheckbox("Use training?", false); // use training data?

        return (gd);

    }

    private double phi(double val1, double val2) { // neighbor interaction,
        // penalize highly for
        // neighbor differences
        return ((val1 - val2) * (val1 - val2) / phi);
    }

    private double psi(double val1, double val2) {
        return ((val1 - val2) * (val1 - val2) / this.psi);
    }

    public void run(ImageProcessor ip) {
        ImageStack GibbsStack = imp.createEmptyStack();
        width = ip.getWidth();
        height = ip.getHeight();
        int slices = imp.getImageStackSize();
        if (slices > 1)
            IJ
                    .error("This only works on single images and not on stacks for now");
        GenericDialog gd = OptionSelector();
        gd.showDialog();
        if (gd.wasCanceled())
            return;

        iters = (int) gd.getNextNumber();

        boolean train = gd.getNextBoolean();

        double[][] X = getMatrix(ip); // observation matrix
        double[][] Y = X.clone();

        if (train) {
            // open a training file
            JFileChooser fc = new JFileChooser();
            fc.showOpenDialog(null);
            trainingfile = fc.getSelectedFile().getPath();
            Opener opener = new Opener();

            ImagePlus trainingImage = opener.openImage(trainingfile);
            ImageProcessor ip2 = trainingImage.getProcessor();
            float[] trainingpixels = (float[]) ip2.getPixels();

            for (int x = 0; x < trainingpixels.length; x++) {
                if (trainingpixels[x] == 1)
                    trainingpixels[x] = 0;
                else
                    trainingpixels[x] = 1;
            }
            ip2.setPixels(trainingpixels);
            trainingImage.show();
            trainingImage.updateAndDraw();
            double[][] Z = getMatrix(ip2); // training image
            learn(Z, X);
        } else {
            phi = gd.getNextNumber();
            theta = gd.getNextNumber();
            psi = gd.getNextNumber();
        }

        FloatProcessor ip3 = MatrixToFloatProcessor(Y);
        GibbsStack.addSlice("", ip3);
        ImagePlus GibbsResult = new ImagePlus("Result", GibbsStack);
        GibbsResult.show();

        for (int iter = 1; iter <= iters; iter++) {
            IJ.showStatus("Gibbs iteration " + iter);

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    double cdf = Math.random();
                    double probs[] = isingprob(Y, x, y, X);
                    if (cdf < probs[0])
                        Y[x][y] = 1;
                    else if (cdf < probs[0] + probs[1])
                        Y[x][y] = 0;
                    else
                        Y[x][y] = -1;
                }
            }
            ip3 = null;
            ip3 = MatrixToFloatProcessor(Y);
            GibbsResult.setProcessor("Iter " + iter + " " + this.theta, ip3);
            GibbsResult.updateAndDraw();
        }
        return;
    }

    public int setup(String arg, ImagePlus imp) {

        this.imp = imp;
        return DOES_32;
    }

    private double theta(double val) {
        return (Math.abs(val) / theta);
    }

}