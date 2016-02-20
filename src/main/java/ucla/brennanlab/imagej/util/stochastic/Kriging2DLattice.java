package ucla.brennanlab.imagej.util.stochastic;

import org.ujmp.core.DenseMatrix2D;
import org.ujmp.core.Matrix;
import org.ujmp.core.SparseMatrix;
import org.ujmp.core.doublematrix.DenseDoubleMatrix;
import org.ujmp.core.doublematrix.SparseDoubleMatrix;
import org.ujmp.core.util.MathUtil;

import java.util.ArrayList;

/**
 * This is a collection of methods for performing Bayesian universal kriging
 * over a 2D lattice, with GMRF. This is kind of a dumb library. Please perform
 * checks outside of this library.
 *
 * @author Josh Chang
 */


public class Kriging2DLattice {
    double a = 1;
    double b = 2;
    int width;
    int height;
    double priorMean, priorPrecision;
    ArrayList<int[]> points;
    ArrayList<Double> measurements;

    int correlationLength = 2;

    public Kriging2DLattice(double priorMean, double priorPrecision, int width, int height){
        this.priorMean = priorMean;
        this.priorPrecision = priorPrecision;
        this.width = width;
        this.height = height;
        this.points = new ArrayList<int[]>();
        this.measurements = new ArrayList<Double>();
    }

    public ArrayList<int[]> getPoints(){
        return this.points;
    }

    public ArrayList<Double> getMeasurements(){
        return this.measurements;
    }

    /**
     * Perform inference using the points in mask. Returns the values
     * and the covariance matrix
     * @param mask   mask for which of the points we want to use
     * @return some matrices: fitted (\hat\beta Eq 6), A_\beta from Eq 7
     */
    public ArrayList<double[][]> infer(boolean[] mask){
        assert(mask.length == this.points.size());
        ArrayList<int[]> points = new ArrayList<int[]>();
        ArrayList<Double> measurements = new ArrayList<Double>();
        int npoints = 0;
        for(int i = 0; i<this.points.size();i++){
            if(mask[i]){
                npoints++;
                points.add(this.points.get(i));
                measurements.add(this.measurements.get(i));
            }
        }
        // populate the precision matrix
        Matrix R0 = SparseDoubleMatrix.Factory.zeros(npoints,npoints);
        double R0sum = 0; // 1^t R0 1
        double[] R0colsums = new double[npoints];
        double tmp;
        double V0R0V0 = 0;

        for (int column = 0; column < npoints; column++) {
            R0colsums[column] = 0;
            for (int row = 0; row < npoints; row++) {

                int x1 = points.get(row)[0];
                int y1 = points.get(row)[1];
                int x2 = points.get(column)[0];
                int y2 = points.get(column)[1];
                if(Math.abs(x1-x2)<=2 && Math.abs(y1-y2)<=1) {
                    tmp = rueGaussianPrecision(x1, y1,x2, y2);
                    R0.setAsDouble(tmp, row, column);
                    R0sum += tmp;
                    R0colsums[column] += tmp;
                    V0R0V0 += tmp * measurements.get(row) * measurements.get(column);
                }
            }
        }
        double[][] betahat = new double[][]{{0.0}};
        double[][] residuals = new double[npoints][1];
        double betahatprefactor = Math.pow(R0sum+this.priorPrecision,-1);

        double betahatsum2 = 0;
        double tmp2;
        for(int i=0; i<npoints; i++ ){
            tmp2 = betahatprefactor*(R0colsums[i]*measurements.get(i)+this.priorPrecision*this.priorMean);
            betahat[0][0] += tmp2;
            betahatsum2+= tmp2*tmp2;
            residuals[i][0] = measurements.get(i) - betahat[0][0];
        }


        double[][] Abeta = new double[1][1];
        Abeta[0][0] = betahatprefactor/(this.a+npoints)*
                (2*this.b+this.priorPrecision*this.priorMean*this.priorMean + V0R0V0 - betahatsum2/betahatprefactor
                );

        ArrayList<double[][]> output = new ArrayList<double[][]>();
        output.add(betahat);
        output.add(Abeta);
        output.add(residuals);

        double arraypoints[][] = new double[npoints][2];
        for(int i=0;i<npoints;i++) {
            arraypoints[i][0] = points.get(i)[0];
            arraypoints[i][1] = points.get(i)[1];
        }

        output.add(arraypoints);
        return output;
    }

    /**
     * Infer based on all of the current points available
     * @return  ArrayList of estimate, variance, points
     */
    public ArrayList<double[][]> infer(){
        boolean[] mask = new boolean[this.measurements.size()];
        for(int i = 0; i<this.measurements.size();i++){
            mask[i] = true;
        }
        return this.infer(mask);
    }

    /**
     * Solve for Eq 8, 9 - the prediction and the variance matrix
     * @param inferred   The inferred values from calling this.infer()
     * @param newpoints  The new points at which to predict based on the inferred values
     * @return Eq 8, 9
     */
    public ArrayList<double[][]> predict(ArrayList<double[][]> inferred, ArrayList<int[]> newpoints){
        ArrayList<double[][]> prediction = new ArrayList<double[][]>();

        double[][] betahat = inferred.get(0);
        double[][] Abeta = inferred.get(1);
        double[][] residuals = inferred.get(2);
        double[][] originalpoints = inferred.get(3);

        double[][] Vhat = new double[newpoints.size()][1];

        int n = originalpoints.length;
        int m = newpoints.size();

        double[][] newpointsarray = new double[newpoints.size()][2];
        for(int i=0;i<m;i++){
            newpointsarray[i][0] = newpoints.get(i)[0];
            newpointsarray[i][1] = newpoints.get(i)[1];
        }

        SparseMatrix U = makePrecisionMatrix(newpointsarray,originalpoints);
        SparseMatrix R = makePrecisionMatrix(newpointsarray);
        SparseMatrix R0 = makePrecisionMatrix(originalpoints);

        Matrix residualMatrix = DenseDoubleMatrix.Factory.zeros(n,1);
        for(int i=0;i<n;i++){
            residualMatrix.setAsDouble(residuals[i][0],i,0);
        }
        Matrix perturbation = R.solve(U.mtimes(residualMatrix));

        for(int i=0;i<m;i++){
            Vhat[i][0] = betahat[0][0] - perturbation.getAsDouble(i,0); // the mean
        }

        Matrix X0 = DenseDoubleMatrix.Factory.zeros(n,1).plus(1);
        Matrix X = DenseDoubleMatrix.Factory.zeros(m,1).plus(1);

        Matrix tmp1 = R.solve(U.mtimes(X0)).plus(X);

        Matrix adjustment = X0.transpose().mtimes(R0).mtimes(X0).plus(this.priorPrecision);

        Matrix Amatrix =
                tmp1.times(Abeta[0][0]).mtimes(tmp1.transpose())
                        .plus(R.inv().times(adjustment.getAsDouble(0,0)).times(Abeta[0][0]));

        // Amatrix.showGUI();
        double[][] A = Amatrix.toDoubleArray();

        prediction.add(Vhat);
        prediction.add(A);
        return prediction;
    }

    public ArrayList<double[][]> sample(ArrayList<double[][]> prediction,int numsamples){
        double[][] estimate = prediction.get(0);
        double[][] A = prediction.get(1);
        long Npoints = estimate.length;

        Matrix P = DenseMatrix2D.Factory.zeros(estimate.length,estimate.length);
        for(long i=0;i<Npoints;i++){
            for(long j=0;j<Npoints;j++){
                P.setAsDouble(A[(int)i][(int)j],i,j);
            }
        }
        Matrix L = P.chol(); // estimate.length x estimate.length

        // create a random Gaussian matrix of white noise
        Matrix Z = DenseMatrix2D.Factory.zeros(estimate.length,numsamples); // estimate.length x numsamples
        for(int j=0;j<numsamples;j++){
            for(int i=0;i<estimate.length;i++){
                Z.setAsDouble(MathUtil.nextGaussian(),i,j);
            }
        }
        Matrix LZ = L.mtimes(Z); // estimate.length x numsamples

        ArrayList<double[][]> samples = new ArrayList<double[][]>(numsamples);

        for(int i=0;i<numsamples;i++){
            double[][] sample = new double[estimate.length][1];
            for(int j = 0; j<estimate.length;j++){
                sample[j][0] = Math.max(estimate[j][0] + LZ.getAsDouble(j,i),0.01);
            }
            samples.add(sample);
        }
        return samples;
    }

    /**
     * From Table 2 in Rue and Tjelmeland
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    private double rueGaussianPrecision(int x1, int y1, int x2, int y2) {
        double precision = Math.pow(0.262, -2);
        if (x1 == x2) {
            if (y1 == y2)
                return precision;
            else if (Math.abs(y1 - y2) == 1)
                return -0.287 * precision;
            else if (Math.abs(y1 - y2) == 2)
                return 0.220 * precision;
            else
                return 0;
        } else if (Math.abs(x1 - x2) == 1) {
            if (y1 == y2)
                return -0.287 * precision;
            else if (Math.abs(y1 - y2) == 1)
                return -0.013 * precision;
            else if (Math.abs(y1 - y2) == 2)
                return -0.128 * precision;
            else
                return 0;
        } else if (Math.abs(x1 - x2) == 2) {
            if (y1 == y2)
                return 0.220 * precision;
            else if (Math.abs(y1 - y2) == 1)
                return -0.128 * precision;
            else if (Math.abs(y1 - y2) == 2)
                return 0.087 * precision;
            else
                return 0;
        } else
            return 0;
    }


    private SparseMatrix makePrecisionMatrix(double[][] points1, double points2[][]){
        int n = points1.length;
        int m = points2.length;
        SparseMatrix out = SparseDoubleMatrix.Factory.zeros(n,m);
        double x1, y1, x2, y2;
        for(int j = 0;j<m; j++){
            for(int i=0;i<n; i++){
                x1 = points1[i][0];
                y1 = points1[i][1];
                x2 = points2[j][0];
                y2 = points2[j][1];
                if(Math.abs(x1-x2)<=2 && Math.abs(y2-y1)<=2){
                    out.setAsDouble(rueGaussianPrecision((int) x1, (int) y1, (int) x2, (int) y2),i,j);
                }
            }
        }
        return out;
    }

    private SparseMatrix makePrecisionMatrix(double[][] points1){
        return makePrecisionMatrix(points1,points1);
    }

    public void addObservation(int x, int y, double value){
        this.points.add(new int[]{x,y});
        this.measurements.add(value);
    }

}
