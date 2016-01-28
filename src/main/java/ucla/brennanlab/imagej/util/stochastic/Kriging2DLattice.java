package ucla.brennanlab.imagej.util.stochastic;

import org.ujmp.core.DenseMatrix2D;
import org.ujmp.core.Matrix;
import org.ujmp.core.doublematrix.SparseDoubleMatrix;

import java.util.ArrayList;

/**
 * This is a collection of methods for performing Bayesian universal kriging
 * over a 2D lattice, with GMRF,
 *
 * @author Josh Chang
 */


public class Kriging2DLattice {
    double a = 1;
    double b = 2;
    long n = 0; // number of inputted points


    private int correlationLength = 5;
    /**
     * When calculating beta for large lattice
     */

    private Matrix V0; // n x 1 matrix
    private ArrayList<int[]> locations; // nx2 matrix, should actually be integer but
    // whatever
    private Matrix X0; // n x p matrix
    // F - Xbeta
    private Matrix residuals;
    private Matrix beta0;
    private Matrix betaHat;
    private Matrix betaHatCovariance;
    /**
     * Some intermediate computation that we will use
     */

    private Matrix X0R0X0;

    // X_0^tR_0X_0+V^{-1} pxp
    private Matrix X0R0X0plusVinv;

    // (X_0^tR_0X_0+V^{-1})^{-1}
    private Matrix invX0R0X0plusVinv;

    private Matrix R0;
    private Matrix P; // P
    //	private SparseMatrix withinPrecision;
    private Matrix predicted;
    private Matrix predictionVariance;
    private double bw;
    private Matrix predictionRootVariance; // cholesky decomposition square root

    /**
     *
     * @param beta0
     * @param P
     */
    public Kriging2DLattice(double beta0, double P){
        this.P = DenseMatrix2D.Factory.zeros(1,1);
        this.P.setAsDouble(P,0,0);

        this.beta0 = DenseMatrix2D.Factory.zeros(1,1);
        this.beta0.setAsDouble(beta0,0,0);
    }

    /**
     *
     * @param beta0
     * @param P
     */
    public Kriging2DLattice(Matrix beta0, Matrix P){
        if (P.getRowCount() != beta0.getRowCount())
            throw new IllegalArgumentException(
                    "Matrix2D dimensions must be compatible");
        this.P = P;
        this.beta0 = beta0;
    }

    public void addObservations(int[][] locations, double[][] X, double[] V){
        ArrayList<int[] > loc = new ArrayList<int[]>(locations.length);
        Matrix X1 = DenseMatrix2D.Factory.zeros(X.length,X[0].length);
        Matrix V1 = DenseMatrix2D.Factory.zeros(V.length,1);

        for(int i=0;i<locations.length; i++){
            loc.add(new int[]{locations[i][0], locations[i][1]});
            V1.setAsDouble(V[i],i,0);
            for(int j=0;j<X[0].length;j++){
                X1.setAsDouble(X[i][j],i,j);
            }
        }

        for(int i=0; i<V.length; i++)

        this.addObservations(loc,X1,V1);
    }

    /**
     * Add new observations and update the Kriging model using block matrix
     * linear algebra
     *
     * @param locations  ArrayList of Point2D
     * @param X nxp matrix of X0
     * @param V   nx1 matrix of V0
     */
    public void addObservations(ArrayList<int[]> locations, Matrix X,
                                Matrix V) {

        if(n==0){
            if(X.getColumnCount() != beta0.getColumnCount()){
                throw new IllegalArgumentException(
                        "Matrix2D dimensions must be compatible");
            }
            if (locations.size() != V.getRowCount())
                throw new IllegalArgumentException(
                        "Matrix2D dimensions for locations and V0 incompatible");
            if (locations.size() != X.getRowCount())
                throw new IllegalArgumentException(
                        "Matrix2D dimensions for X0 and locations incompatible");

            this.X0 = X;
            this.locations = locations;
            this.V0 = V;
            this.n = V.getRowCount();

            /**
             * 5x5 GMRF
             */
            this.R0 = SparseDoubleMatrix.Factory
                    .zeros(this.n, this.n);
            // Set the correlation matrix
            for (int row = 0; row < n; row++) {
                for (int column = 0; column < n; column++) {
                    R0.setAsDouble(rueGaussianPrecision(locations.get(row)[0], locations.get(row)[1],
                            locations.get(column)[0], locations.get(
                                    column)[1]), row, column);
                }
            }
        }else{
            int nobs1 = (int) this.V0.getRowCount();
            int nobs2 = (int) V.getRowCount();
            int nobs =  nobs1 + nobs2;
            Matrix newresponses = DenseMatrix2D.Factory.zeros(nobs,1);
            int ncovariates = (int) this.X0.getColumnCount();
            Matrix newcovariates = DenseMatrix2D.Factory.zeros(nobs,ncovariates);

            for(int j=0;j<nobs1;j++){
                newresponses.setAsDouble(this.V0.getAsDouble(j,0),j,0);
                for(int k=0;k<ncovariates;k++){
                    newcovariates.setAsDouble(this.X0.getAsDouble(j,k),j,k);
                }
            }
            for(int j=0;j<nobs2;j++){
                newresponses.setAsDouble(V.getAsDouble(j,0),j+nobs1,0);
                for(int k=0;k<ncovariates;k++){
                    newcovariates.setAsDouble(X.getAsDouble(j,k),j+nobs1,k);
                }
            }

            this.V0 = V;
            this.X0 = newcovariates;
            Matrix withinPrecision = SparseDoubleMatrix.Factory.zeros(locations
                    .size(), locations.size());

            for(int j=0; j<locations.size();j++){
                this.locations.add(locations.get(j));
            }

            for (int column = 0; column < nobs; column++) {
                for (int row = 0; row < nobs; row++) {
                    withinPrecision.setAsDouble(rueGaussianPrecision(this.locations
                                    .get(row)[0], this.locations.get(row)[1],
                            this.locations.get(column)[0], this.locations.get(
                                    column)[1]), row, column);
                }
            }
            this.R0 = withinPrecision;
        }




        computeBeta();

    }

    public void clean() {
        this.X0 = null;
        this.locations = null;
        this.R0 = null; // R_0
    }

    /**
     * Compute beta and its covariance matrix
     */
    public void computeBeta() {

        X0R0X0plusVinv = quadraticInnerNorm(R0,
                X0).plus(P);
        invX0R0X0plusVinv = X0R0X0plusVinv.inv();

        betaHat = invX0R0X0plusVinv.mtimes(X0.transpose()
                .mtimes(R0).mtimes(V0).plus(
                        P.mtimes(beta0)));

        bw = (2 * b + V0
                .transpose()
                .mtimes(R0)
                .mtimes(V0)
                .minus(
                        quadraticInnerNorm(X0R0X0plusVinv, betaHat)
                                .plus(
                                        quadraticInnerNorm(this.P,
                                                this.beta0))).getAsDouble(
                        0, 0));

        betaHatCovariance = invX0R0X0plusVinv.times(bw / (n + 2 * a));

        residuals = V0.minus(this.X0.mtimes(this.betaHat));
    }



    /**
     * @return
     */
    public Matrix getBeta() {
        return this.betaHat;
    }

    /**
     * @return
     */
    public Matrix getBetaVariance() {
        return this.betaHatCovariance;
    }

    public Matrix getPredicted() {
        return predicted;
    }

    public void setPredicted(Matrix predicted) {
        this.predicted = predicted;
    }

    public Matrix getPredictionRootVariance() {
        return predictionRootVariance;
    }

    public void setPredictionRootVariance(Matrix predictionRootVariance) {
        this.predictionRootVariance = predictionRootVariance;
    }

    /**
     * Interpolate the field at these points, and give covariance matrix while
     * at it (Normal approximation)
     *
     * @param locations
     */
    public void interpolate(ArrayList<int[]> locations, Matrix covariates) {

        assert (locations.size() == covariates.getColumnCount());

        /**
         * First reshape the observations to cull out far-away points. This step
         * improves efficiency, and also improves numerical stability because we
         * remove a lot of (close to zero) eigenvalues that can screw up a
         * Cholesky decomposition. We do this in a pretty stupid and lazy way
         * however using signed distance functions, just out of laziness
         */

        // Compare this.locations and locations
		/*
		 * Matrix culledLocations; Matrix culledResiduals; Matrix
		 * culledCovariates;
		 *
		 * // Maximum X and Y in lattice double maxX, maxY; // Minimum X and Y
		 * in a lattice double minX,minY; int Xrange = (int) (maxX-minX); int
		 * Yrange = (int) (maxY-minY);
		 *
		 *
		 *
		 * for(int i=0; i<this.locations.getRowCount();i++){
		 *
		 * }
		 */

        Matrix partialPrecision = SparseDoubleMatrix.Factory.zeros(locations
                .size(), this.locations.size());
        long columns = this.locations.size();
        long rows = locations.size();
        // Compute cross-correlations
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                partialPrecision.setAsDouble(rueGaussianPrecision(locations
                                .get(row)[0], locations.get(row)[1],
                        this.locations.get(column)[0], this.locations
                                .get(column)[1]), row, column);

            }
        }
        // compute within correlations
        Matrix withinPrecision = SparseDoubleMatrix.Factory.zeros(locations
                .size(), locations.size());

        for (int column = 0; column < rows; column++) {
            for (int row = 0; row < rows; row++) {
                withinPrecision.setAsDouble(rueGaussianPrecision(locations
                                .get(row)[0], locations.get(row)[1],
                        locations.get(column)[0], locations.get(
                                column)[1]), row, column);
            }
        }

        // Use linear algebra solver
        Matrix newResiduals = withinPrecision.solve(
                partialPrecision.mtimes(residuals)).times(-1);
        this.setPredicted(newResiduals.plus(covariates.mtimes(this.betaHat)));

        this.predictionVariance = quadraticOuterNorm(
                this.betaHatCovariance,
                covariates.minus(withinPrecision.solve(partialPrecision
                        .mtimes(this.X0)))).plus(
                withinPrecision.times(bw / (n + 2 * this.a)));

        this.setPredictionRootVariance(predictionVariance.chol());
        //
    }

    /**
     * Return B^t A B
     *
     * @param A square matrix of dimension nxn
     * @param B matrix of dimension nxm
     * @return B^tAB
     */
    public Matrix quadraticInnerNorm(Matrix A, Matrix B) {
        if (A.getRowCount() != A.getColumnCount())
            throw new IllegalArgumentException("A must be square");
        if (A.getRowCount() != B.getRowCount())
            throw new IllegalArgumentException("B and A are not compatible");
        return B.transpose().mtimes(A.mtimes(B));
    }

    /**
     * Return B A B^t
     *
     * @param A
     * @param B
     * @return BAB^t
     */
    public Matrix quadraticOuterNorm(Matrix A, Matrix B) {
        if (A.getRowCount() != A.getColumnCount())
            throw new IllegalArgumentException("A must be square");
        if (A.getRowCount() != B.getColumnCount())
            throw new IllegalArgumentException("B and A are not compatible");
        return B.mtimes(A).mtimes(B.transpose());
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

    /**
     * Subtract two matrices
     *
     * @param A
     * @param B
     * @return The difference A-B
     */
    public Matrix subtractMatrices(Matrix A, Matrix B) {
        if (A.getRowCount() != B.getRowCount())
            throw new IllegalArgumentException("Matrix2D dimensions must agree");
        if (A.getColumnCount() != B.getColumnCount())
            throw new IllegalArgumentException("Matrix2D dimensions must agree");

        return A.minus(B);
    }




}
