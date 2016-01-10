package ucla.brennanlab.imagej.util.stochastic;

import org.ujmp.core.Matrix;
import org.ujmp.core.doublematrix.SparseDoubleMatrix;
import ucla.brennanlab.imagej.util.levelsets.ImplicitShape2D;


/**
 * This is a collection of methods for performing Bayesian universal kriging
 * over a 2D lattice, with GMRF,
 *
 * @author Josh Chang
 */

public class Kriging2DLattice {
    double a = 1;
    double b = 2;
    long n; // dimensions of the kriging signedDistance
    private Matrix response; // n x 1 matrix
    private Matrix locations; // nx2 matrix, should actually be integer but
    // whatever
    private Matrix covariates; // n x p matrix

    private Matrix betaPrior;
    private Matrix betaHat;
    private Matrix betaHatCovariance;

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

    public Kriging2DLattice(Matrix covariates, Matrix locations,
                            Matrix response, Matrix betaprior, Matrix betaPrecision) {
        if (locations.getRowCount() != response.getRowCount())
            throw new IllegalArgumentException(
                    "Matrix2D dimensions for locations and response incompatible");
        if (locations.getRowCount() != covariates.getRowCount())
            throw new IllegalArgumentException(
                    "Matrix2D dimensions for covariates and locations incompatible");
        if (betaPrecision.getRowCount() != betaprior.getRowCount()
                || covariates.getColumnCount() != betaprior.getColumnCount())
            throw new IllegalArgumentException(
                    "Matrix2D dimensions must be compatible");
        this.betaPrecision = betaPrecision;
        this.covariates = covariates;
        this.locations = locations;
        this.response = response;
        this.n = response.getRowCount();

        /**
         * 5x5 GMRF
         */
        this.spatialPrecision = SparseDoubleMatrix.Factory
                .zeros(this.n, this.n);
        this.betaPrior = betaprior;
        // Set the correlation matrix
        for (long row = 0; row < n; row++) {
            for (long column = 0; column < n; column++) {
                spatialPrecision.setAsDouble(rueGaussianPrecision(locations
                                .getAsInt(row, 0), locations.getAsInt(row, 1),
                        locations.getAsInt(column, 0), locations.getAsInt(
                                column, 1)), row, column);
            }
        }
        // spatialPrecision.showGUI();

    }

    /**
     * Approximately!!!
     * Add observations to estimate a new betaHat. We do this in a way
     * that is completely wrong, by just assuming independence and piecing together
     * the results to form the new estimate, but we can fix this later TODO, fix this, but
     * this won't really change the solution much
     *
     * @param locations  nx2 matrix of locations
     * @param covariates nxp matrix of covariates
     * @param response   nx1 matrix of response
     */
    public void addObservations(Matrix locations, Matrix covariates,
                                Matrix response) {


        Matrix withinPrecision = SparseDoubleMatrix.Factory.zeros(locations
                .getRowCount(), locations.getRowCount());

        long rows = locations.getRowCount();

        for (int column = 0; column < rows; column++) {
            for (int row = 0; row < rows; row++) {
                withinPrecision.setAsDouble(rueGaussianPrecision(locations
                                .getAsInt(row, 0), locations.getAsInt(row, 1),
                        locations.getAsInt(column, 0), locations.getAsInt(
                                column, 1)), row, column);
            }
        }
/*		
        Matrix partialunscaledCondBetaVarianceInv = quadraticInnerNorm(withinPrecision,
				covariates).plus(betaPrecision);

		Matrix partialunscaledCondBetaVariance = partialunscaledCondBetaVarianceInv.inv();
		Matrix partialbetaHat = partialunscaledCondBetaVariance.mtimes(covariates.transpose()
				.mtimes(withinPrecision).mtimes(response).plus(
						betaPrecision.mtimes(betaPrior)));
		
		double partialbw = (2 * b + response
				.transpose()
				.mtimes(withinPrecision)
				.mtimes(response)
				.minus(
						quadraticInnerNorm(partialunscaledCondBetaVarianceInv, partialbetaHat)
								.plus(
										quadraticInnerNorm(this.betaPrecision,
												this.betaPrior))).getAsDouble(
						0, 0));
*/
        //	Matrix partialbetaHatCovariance = partialunscaledCondBetaVariance.times(partialbw / (n + 2 * a));

        //	Matrix partialresiduals = response.minus(covariates.mtimes(partialbetaHat));

        // Now augment the original estimates

    }

    public void clean() {
        this.covariates = null;
        this.locations = null;
        this.spatialPrecision = null;
    }

    /**
     * Compute beta and its covariance matrix
     */
    public void computeBeta() {

        unscaledCondBetaVarianceInv = quadraticInnerNorm(spatialPrecision,
                covariates).plus(betaPrecision);
        unscaledCondBetaVariance = unscaledCondBetaVarianceInv.inv();

        betaHat = unscaledCondBetaVariance.mtimes(covariates.transpose()
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

        residuals = response.minus(this.covariates.mtimes(this.betaHat));
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
                .getRowCount(), this.locations.getRowCount());
        long columns = this.locations.getRowCount();
        long rows = locations.getRowCount();
        // Compute cross-correlations
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                partialPrecision.setAsDouble(rueGaussianPrecision(locations
                                .getAsInt(row, 0), locations.getAsInt(row, 1),
                        this.locations.getAsInt(column, 0), this.locations
                                .getAsInt(column, 1)), row, column);

            }
        }
        //partialPrecision.showGUI();
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
                        .mtimes(this.covariates)))).plus(
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
