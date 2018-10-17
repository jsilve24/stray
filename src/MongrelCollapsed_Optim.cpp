#include <MatrixAlgebra.h>
#include <MongrelCollapsed.h>
#include <AdamOptim.h>
#include <AdamOptimPerturb.h> // optional not fully implemented yet (or helpful)
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;

//' Function to Optimize the Collapsed Mongrel Model
//' 
//' See details for model. Should likely be followed by function 
//' \code{\link{uncollapseMongrelCollapsed}}. Notation: \code{N} is number of samples,
//' \code{D} is number of multinomial categories, and \code{Q} is number
//' of covariates. 
//' 
//' @param Y D x N matrix of counts
//' @param upsilon (must be > D)
//' @param ThetaX D-1 x N matrix formed by Theta*X (Theta is Prior mean 
//'    for regression coefficients) 
//' @param K D-1 x D-1 precision matrix (inverse of Xi)
//' @param A N x N precision matrix given by (I_N + X*Gamma*X')^{-1}]
//' @param init D-1 x N matrix of initial guess for eta used for optimization
//' @param n_samples number of samples for Laplace Approximation (=0 very fast
//'    as no inversion or decomposition of Hessian is required)
//' @param calcGradHess if n_samples=0 should Gradient and Hessian 
//'   still be calculated using closed form solutions?
//' @param b1 (ADAM) 1st moment decay parameter (recommend 0.9) "aka momentum"
//' @param b2 (ADAM) 2nd moment decay parameter (recommend 0.99 or 0.999)
//' @param step_size (ADAM) step size for descent (recommend 0.001-0.003)
//' @param epsilon (ADAM) parameter to avoid divide by zero
//' @param eps_f (ADAM) normalized function improvement stopping criteria 
//' @param eps_g (ADAM) normalized gradient magnitude stopping criteria
//' @param max_iter (ADAM) maximum number of iterations before stopping
//' @param verbose (ADAM) if true will print stats for stopping criteria and 
//'   iteration number
//' @param verbose_rate (ADAM) rate to print verbose stats to screen
//' @param decomp_method decomposition of hessian for Laplace approximation
//'   'eigen' (more stable, slower, default) or 'cholesky' (less stable, faster)
//' @param eigvalthresh threshold for negative eigenvalues in 
//'   decomposition of negative inverse hessian (should be <=0)
//' @param no_error if true will throw hessian warning rather than error if 
//'   not positive definite. 
//' @param jitter (default: 0) if >0 then adds that factor to diagonal of Hessian 
//' before decomposition (to improve matrix conditioning)
//'   
//' @details Notation: Let Z_j denote the J-th row of a matrix Z.
//' Model:
//'    \deqn{Y_j ~ Multinomial(Pi_j)}
//'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
//'    \deqn{Eta ~ T_{D-1, N}(upsilon, Theta*X, K^{-1}, A^{-1})}
//' Where A = (I_N + X * Gamma * X')^{-1}, K^{-1} = Xi is a (D-1)x(D-1) covariance 
//' matrix, Gamma is a Q x Q covariance matrix, and Phi^{-1} is ALRInv_D 
//' transform. 
//' 
//' Gradient and Hessian calculations are fast as they are computed using closed
//' form solutions. That said, the Hessian matrix can be quite large 
//' \[N*(D-1) x N*(D-1)\] and storage may be an issue. 
//' 
//' Note: Warnings about large negative eigenvalues can either signal 
//' that the optimizer did not reach an optima or (more commonly in my experience)
//' that the prior / degrees of freedom for the covariance (given by parameters
//' \code{upsilon} and \code{K}) were too specific and at odds with the observed data.
//' If you get this warning try the following. 
//' 1. Try restarting the optimization using a different initial guess for eta
//' 2. Try decreasing (or even increasing )\code{step_size} (by increments of 0.001 or 0.002) 
//'   and increasing \code{max_iter} parameters in optimizer. Also can try 
//'   increasing \code{b1} to 0.99 and decreasing \code{eps_f} by a few orders
//'   of magnitude
//' 3. Try relaxing prior assumptions regarding covariance matrix. (e.g., may want
//' to consider decreasing parameter \code{upsilon} closer to a minimum value of 
//' D)
//' 4. Try adding small amount of jitter (e.g., set \code{jitter=1e-5}) to address
//'   potential floating point errors. 
//' @return List containing (all with respect to found optima)
//' 1. LogLik - Log Likelihood of collapsed model (up to proportionality constant)
//' 2. Gradient - (if \code{calcGradHess}=true)
//' 3. Hessian - (if \code{calcGradHess}=true)
//' 4. Pars - Parameter value of eta at optima
//' 5. Samples - (D-1) x N x n_samples array containing posterior samples of eta 
//'   based on Laplace approximation (if n_samples>0)
//' @md 
//' @export
//' @name optimMongrelCollapsed
//' @references S. Ruder (2016) \emph{An overview of gradient descent 
//' optimization algorithms}. arXiv 1609.04747
//' @seealso \code{\link{uncollapseMongrelCollapsed}}
//' @examples
//' sim <- mongrel_sim()
//' 
//' # Fit model for eta
//' fit <- optimMongrelCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$K, 
//'                              sim$A, random_mongrel_init(sim$Y))  
// [[Rcpp::export]]
List optimMongrelCollapsed(const Eigen::ArrayXXd Y, 
               const double upsilon, 
               const Eigen::MatrixXd ThetaX, 
               const Eigen::MatrixXd K, 
               const Eigen::MatrixXd A, 
               Eigen::MatrixXd init, 
               int n_samples=2000, 
               bool calcGradHess = true,
               double b1 = 0.9,         
               double b2 = 0.99,        
               double step_size = 0.003, // was called eta in ADAM code
               double epsilon = 10e-7, 
               double eps_f=1e-10,       
               double eps_g=1e-4,       
               int max_iter=10000,      
               bool verbose=false,      
               int verbose_rate=10,
               String decomp_method="eigen",
               double eigvalthresh=0, 
               bool no_error=false, 
               double jitter=0){  
  int N = Y.cols();
  int D = Y.rows();
  MongrelCollapsed cm(Y, upsilon, ThetaX, K, A);
  Map<VectorXd> eta(init.data(), init.size()); // will rewrite by optim
  double nllopt; // NEGATIVE LogLik at optim
  List out(5);
  out.names() = CharacterVector::create("LogLik", "Gradient", "Hessian",
            "Pars", "Samples");
  
  // Pick optimizer (ADAM - without perturbation appears to be best)
  //   ADAM with perturbations not fully implemented
  //int status = Numer::optim_lbfgs(cm, eta, nllopt);
  int status = adam::optim_adam(cm, eta, nllopt, b1, b2, step_size, epsilon, 
                                eps_f, eps_g, max_iter, verbose, verbose_rate); 
  //int status = adamperturb::optim_adam(cm, eta, nllopt); 

  if (status<0)
    Rcpp::warning("Max Iterations Hit, May not be at optima");
  Map<MatrixXd> etamat(eta.data(), D-1, N);
  out[0] = -nllopt; // Return (positive) LogLik
  out[3] = etamat;
  
  if (n_samples > 0 || calcGradHess){
    if (verbose) Rcout << "Allocating for Hessian" << std::endl;
    MatrixXd hess(N*(D-1), N*(D-1));
    VectorXd grad(N*(D-1));
    if (verbose) Rcout << "Calculating Hessian" << std::endl;
    grad = cm.calcGrad(); // should have eta at optima already
    hess = cm.calcHess(); // should have eta at optima already
    if (jitter > 0)       // potentially add jitter to improve conditioning
      hess.diagonal().array() += -jitter;
    if (verbose) Rcout << "Copying Hessian" << std::endl;
    out[1] = grad;
    out[2] = hess;
    
    if (n_samples>0){
      // Laplace Approximation
      if (decomp_method == "eigen"){
        // Eigenvalue Decomposition 
        if (verbose) Rcout << "Decomposing Hessian" << std::endl;
        Eigen::SelfAdjointEigenSolver<MatrixXd> eh(-hess); // negative hessian
        if (verbose) Rcout << "Inverting Eigenvalues" << std::endl;
        VectorXd evalinv(eh.eigenvalues().array().inverse().matrix());
        int excess=0;
        for (int i=1; i<N*(D-1); i++){
          if (evalinv(i) < eigvalthresh) {
            excess++;
          }
        }
        if (excess > 0){
          Rcpp::warning("Some eigenvalues are below minimum threshold");
          Rcout << "Eigenvalues" << evalinv.transpose() << std::endl;
          return out;
        }
        int pos = 0;
        for (int i = N*(D-1)-1; i>=0; i--){
          if (evalinv(pos) > 0)
            pos++;
        }
        if (pos < N*(D-1)-1) {
          Rcpp::warning("Some small negative eigenvalues are being chopped");
          Rcout << N*(D-1)-pos << " out of " << N*(D-1) <<
            " passed eigenvalue threshold"<< std::endl;
        }
        MatrixXd hesssqrt(N*(D-1), pos);
        if (verbose) Rcout << "Calculating Square Root" << std::endl;
        hesssqrt = eh.eigenvectors().rightCols(pos)*
          evalinv.tail(pos).cwiseSqrt().asDiagonal(); //V*D^{-1/2}
        // now generate random numbers...
        if (verbose) Rcout << "Sampling Random Numbers" << std::endl;
        NumericVector r(n_samples*pos);
        r = rnorm(n_samples*pos, 0, 1); // using vectorization from Rcpp sugar
        Map<VectorXd> rvec(as<Map<VectorXd> >(r));
        Map<MatrixXd> rmat(rvec.data(), pos, n_samples);
        MatrixXd samp(pos, n_samples);
        samp.noalias() = hesssqrt*rmat;
        samp.colwise() += eta; // add mean of approximation
        IntegerVector d = IntegerVector::create(D-1, N, n_samples);
        NumericVector samples = wrap(samp);
        samples.attr("dim") = d; // convert to 3d array for return to R
        out[4] = samples;
        
      } else if (decomp_method == "cholesky"){
        //Cholesky Decomposition
        Eigen::LLT<MatrixXd> hesssqrt;
        hesssqrt.compute(-hess);
          if (hesssqrt.info() == Eigen::NumericalIssue){
          if (no_error){
            Rcpp::warning("Cholesky of Hessian failed with status status Eigen::NumericalIssue");
          } else if (!no_error){
              Rcpp::stop("Cholesky of Hessian failed with status Eigen::NumericalIssue");
          }
        }
        NumericVector r(n_samples*N*(D-1));
        r = rnorm(n_samples*N*(D-1), 0, 1); // using vectorization from Rcpp sugar
        Map<VectorXd> rvec(as<Map<VectorXd> >(r));
        Map<MatrixXd> rmat(rvec.data(), N*(D-1), n_samples);
        MatrixXd samp(N*(D-1), n_samples);
        samp.noalias() = hesssqrt.matrixL().solve(rmat);  // calculate errors of approximation
        samp.colwise() += eta; // add mean of approximation
        IntegerVector d = IntegerVector::create(D-1, N, n_samples);
        NumericVector samples = wrap(samp);
        samples.attr("dim") = d; // convert to 3d array for return to R
        out[4] = samples;
      } //endif decomp_method
    } // endif n_samples || calcGradHess
  } // endif n_samples 
  return out;
}
