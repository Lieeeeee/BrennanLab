package ucla.brennanlab.imagej.util.stochastic;

import org.ujmp.core.DenseMatrix2D;
import org.ujmp.core.Matrix;
import org.ujmp.core.Matrix2D;
import org.ujmp.core.doublematrix.SparseDoubleMatrix;
import ucla.brennanlab.imagej.util.Point2D;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;

import java.util.ArrayList;
import java.util.HashSet;

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
    private Matrix response; // n x 1 matrix
    private ArrayList<Point2D> locations; // nx2 matrix, should actually be integer but
    // whatever
    private Matrix X0; // n x p matrix
    private Matrix R0; // n x n matrix

    private Matrix betaPrior;
    private Matrix betaHat;
    private Matrix betaHatCovariance; // P
    private int xreduction;
    private int yreduction;
    /**
     * When calculating beta for large lattice
     */

    // F - Xbeta
    private Matrix residuals;

    /**
     * Some intermediate computation that we will use
     */

    // X_0^tR_0X_0+V^{-1} pxp
    private Matrix unscaledCondBetaVarianceInv;

    // (X_0^tR_0X_0+V^{-1})^{-1}
    private Matrix unscaledCondBetaVariance;

    private Matrix spatialPrecision;
    private Matrix betaPrecision;
    //	private SparseMatrix withinPrecision;
    private Matrix predicted;
    private Matrix predictionVariance;
    private double bw;

    private Matrix predictionRootVariance; // cholesky decomposition square root

    /**
     *
     * @param betaPrior
     * @param betaPrecision
     * @param xreduction
     * @param yreduction
     */
    public Kriging2DLattice(double betaPrior, double betaPrecision, int xreduction, int yreduction){
        this.betaPrecision = DenseMatrix2D.Factory.zeros(1,1);
        this.betaPrecision.setAsDouble(betaPrecision,0,0);

        this.betaPrior = DenseMatrix2D.Factory.zeros(1,1);
        this.betaPrior.setAsDouble(betaPrior,0,0);
        this.xreduction = xreduction;
        this.yreduction = yreduction;
    }

    /**
     *
     * @param betaPrior
     * @param betaPrecision
     */
    public Kriging2DLattice(double betaPrior, double betaPrecision){
        this( betaPrior,  betaPrecision,1,1);
    }

    /**
     *
     * @param betaPrior
     * @param betaPrecision
     * @param xreduction
     * @param yreduction
     */
    public Kriging2DLattice(Matrix betaPrior, Matrix betaPrecision, int xreduction, int yreduction){
        if (betaPrecision.getRowCount() != betaPrior.getRowCount())
            throw new IllegalArgumentException(
                    "Matrix2D dimensions must be compatible");
        this.betaPrecision = betaPrecision;
        this.betaPrior = betaPrior;
        this.xreduction = xreduction;
        this.yreduction = yreduction;
    }


    /**
     * Convenience method for adding observations to the kriging model without
     * first creating UJMP matrices
     * @param locations
     * @param covariates
     * @param response
     */
    public void addObservations(int[][] locations, double[][] covariates, double[] response){
        ArrayList<Point2D> loc = new ArrayList<Point2D>(locations.length);
        for(int j=0; j<locations.length;j++){
            loc.add(new Point2D(locations[j]));
        }
        Matrix cov = DenseMatrix2D.Factory.zeros(covariates.length,covariates[0].length);
        for(int j=0;j<covariates.length;j++){
            for(int k=0;k<covariates[0].length;k++){
                cov.setAsDouble(covariates[j][k],j,k);
            }
        }
        Matrix resp = DenseMatrix2D.Factory.zeros(response.length,1);
        for(int j=0; j<response.length;j++){
            resp.setAsDouble(response[j],j,0);
        }

        addObservations(loc,cov,resp);
        computeBeta();
    }

    /**
     * Add new observations and update the Kriging model using block matrix
     * linear algebra
     *
     * @param locations  nx2 matrix of locations
     * @param covariates nxp matrix of X0
     * @param response   nx1 matrix of response
     */
    public void addObservations(ArrayList<Point2D> locations, Matrix covariates,
                                Matrix response) {

        if(n==0){
            if(covariates.getColumnCount() != betaPrior.getColumnCount()){
                throw new IllegalArgumentException(
                        "Matrix2D dimensions must be compatible");
            }
            if (locations.size() != response.getRowCount())
                throw new IllegalArgumentException(
                        "Matrix2D dimensions for locations and response incompatible");
            if (locations.size() != covariates.getRowCount())
                throw new IllegalArgumentException(
                        "Matrix2D dimensions for X0 and locations incompatible");

            this.X0 = covariates;
            this.locations = locations;
            this.response = response;
            this.n = response.getRowCount();

            /**
             * 5x5 GMRF
             */
            this.spatialPrecision = SparseDoubleMatrix.Factory
                    .zeros(this.n, this.n);
            // Set the correlation matrix
            for (int row = 0; row < n; row++) {
                for (int column = 0; column < n; column++) {
                    spatialPrecision.setAsDouble(rueGaussianPrecision(locations
                                    .get(row).x, locations.get(row).y,
                            locations.get(column).x, locations.get(
                                    column).y), row, column);
                }
            }
        }else{
            int nobs1 = (int) this.response.getRowCount();
            int nobs2 = (int) response.getRowCount();
            int nobs =  nobs1 + nobs2;
            Matrix newresponses = DenseMatrix2D.Factory.zeros(nobs,1);
            int ncovariates = (int) this.X0.getColumnCount();
            Matrix newcovariates = DenseMatrix2D.Factory.zeros(nobs,ncovariates);

            for(int j=0;j<nobs1;j++){
                newresponses.setAsDouble(this.response.getAsDouble(j,0),j,0);
                for(int k=0;k<ncovariates;k++){
                    newcovariates.setAsDouble(this.X0.getAsDouble());
                }
            }
            for(int j=0;j<nobs2;j++){
                newresponses.setAsDouble(response.getAsDouble(nobs1+j,0),j,0);

            }


        }

        Matrix withinPrecision = SparseDoubleMatrix.Factory.zeros(locations
                .size(), locations.size());

        long rows = locations.size();

        for (int column = 0; column < rows; column++) {
            for (int row = 0; row < rows; row++) {
                withinPrecision.setAsDouble(rueGaussianPrecision(locations
                                .get(row).x, locations.get(row).y,
                        locations.get(column).x, locations.get(
                                column).y), row, column);
            }
        }

    }

    public void clean() {
        this.X0 = null;
        this.locations = null;
        this.spatialPrecision = null;
    }

    /**
     * Compute beta and its covariance matrix
     */
    public void computeBeta() {

        unscaledCondBetaVarianceInv = quadraticInnerNorm(spatialPrecision,
                X0).plus(betaPrecision);
        unscaledCondBetaVariance = unscaledCondBetaVarianceInv.inv();

        betaHat = unscaledCondBetaVariance.mtimes(X0.transpose()
                .mtimes(spatialPrecision).mtimes(response).plus(
                        betaPrecision.mtimes(betaPrior)));
        bw = (2 * b + response
                .transpose()
                .mtimes(spatialPrecision)
                .mtimes(response)
                .minus(
                        quadraticInnerNorm(unscaledCondBetaVarianceInv, betaHat)
                                .plus(
                                        quadraticInnerNorm(this.betaPrecision,
                                                this.betaPrior))).getAsDouble(
                        0, 0));

        betaHatCovariance = unscaledCondBetaVariance.times(bw / (n + 2 * a));

        residuals = response.minus(this.X0.mtimes(this.betaHat));
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
     * Interpolate the values within a narrow band downwind in the level set
     *
     * @param signedDistance
     * @param Covariates
     * @param lowerbound
     * @param upperbound
     */
    public void interpolate(ImplicitShape2D signedDistance,
                            Matrix[] Covariates, float lowerbound, float upperbound) {
        /**
         * Basically take elements from the band and feed them into
         * interpolate()
         */
    }

    /**
     * Interpolate the field at these points, and give covariance matrix while
     * at it (Normal approximation)
     *
     * @param locations
     */
    public void interpolate(Matrix locations, Matrix covariates) {

        assert (locations.getRowCount() == covariates.getColumnCount());

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
                .getRowCount(), this.locations.size());
        long columns = this.locations.size();
        long rows = locations.getRowCount();
        // Compute cross-correlations
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                partialPrecision.setAsDouble(rueGaussianPrecision(locations
                                .getAsInt(row, 0), locations.getAsInt(row, 1),
                        this.locations.get(column).x, this.locations
                                .get(column).y), row, column);

            }
        }
        // compute within correlations
        Matrix withinPrecision = SparseDoubleMatrix.Factory.zeros(locations
                .getRowCount(), locations.getRowCount());

        for (int column = 0; column < rows; column++) {
            for (int row = 0; row < rows; row++) {
                withinPrecision.setAsDouble(rueGaussianPrecision(locations
                                .getAsInt(row, 0), locations.getAsInt(row, 1),
                        locations.getAsInt(column, 0), locations.getAsInt(
                                column, 1)), row, column);
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

    /**
     * Get the x-y coordinates corresponding to (x,y) in the Coarsened grid
     * @param x
     * @param y
     * @return
     */
    public Point2D getCoarseCoords(int x, int y){
        return new Point2D((int) Math.floor(x/this.xreduction),
                (int)Math.floor(y/this.yreduction));
    }

    public int getXcoord(int x){
        return (int) Math.floor(x/this.xreduction);
    }

    public int getYcoord(int y){
        return (int)Math.floor(y/this.yreduction);
    }

    /**
     *
     * @param locations
     * @return ArrayList of Matrices. The first matrix is the mean and the second is the variance
     */
    public ArrayList<Matrix2D> infer(int[][] locations){

        HashSet<Point2D> uniqueGridLocations = new HashSet<Point2D>();
        for(int i=0;i<locations.length;i++){
            uniqueGridLocations.add(getCoarseCoords(locations[i][0],locations[i][1]));
        }

        // Figure out all of the grid points that are within range of these locations for inference


        // Compute \hat\beta for the grid points

        return null;
    }



}
