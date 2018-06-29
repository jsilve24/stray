#include <MatrixAlgebra.h>
#include <MultimMatTCollapsed.h>
#include <AdamOptim.h>
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;

//' define optimization function
//' @export
// [[Rcpp::export]]
List optimMMTC(const Eigen::ArrayXXd Y, 
               const double upsilon, 
               const Eigen::MatrixXd ThetaX, 
               const Eigen::MatrixXd K, 
               const Eigen::MatrixXd A, 
               Eigen::MatrixXd etainit, 
               int iter=2000, 
               double eigvalthresh=-0.001, // ignored currently
               int numexcessthresh=0,      // ignored currently 
               bool calcGradHess = true){
  int N = Y.cols();
  int D = Y.rows();
  MultimMatTCollapsed cm(Y, upsilon, ThetaX, K, A);
  Map<VectorXd> eta(etainit.data(), etainit.size()); // will rewrite by optim
  double nllopt; // NEGATIVE LogLik at optim
  List out(5);
  out.names() = CharacterVector::create("LogLik", "Gradient", "Hessian",
            "Pars", "Samples");
  //int status = Numer::optim_lbfgs(cm, eta, nllopt);
  int status = adam::optim_adam(cm, eta, nllopt);
  if (status<0)
    Rcpp::stop("failed to converge");
  Map<MatrixXd> etamat(eta.data(), D-1, N);
  out[0] = -nllopt; // Return (positive) LogLik
  out[3] = etamat;
  
  if (iter > 0 || calcGradHess){
    MatrixXd hess(N*(D-1), N*(D-1));
    VectorXd grad(N*(D-1));
    grad = cm.calcGrad(); // should have eta at optima already
    hess = cm.calcHess(); // should have eta at optima already
    out[1] = grad;
    out[2] = hess;
    
    if (iter>0){
      // Laplace Approximation
      Eigen::SelfAdjointEigenSolver<MatrixXd> eh(-hess); // negative hessian
      VectorXd evalinv(eh.eigenvalues().array().inverse().matrix());
      int excess=0;
      for (int i=1; i<N*(D-1); i++){
        if (evalinv(i) < eigvalthresh) {
          excess++;
        }
      }
      if (excess > numexcessthresh){
        Rcout << evalinv.head(10).transpose() << std::endl;
        Rcpp::stop("To many eigenvalues are below minimum threshold");
      }
      int pos = 0;
      for (int i = N*(D-1)-1; i>=0; i--){
        if (evalinv(pos) > 0)
          pos++;
      }
      if (pos < N*(D-1)) {
        warning("Some small negative eigenvalues are being chopped");
        Rcout << N*(D-1)-pos << " out of " << N*(D-1) <<
          " passed eigenvalue threshold"<< std::endl;
      }
      MatrixXd hesssqrt(N*(D-1), pos);
      hesssqrt = eh.eigenvectors().rightCols(pos)*
        evalinv.tail(pos).cwiseSqrt().asDiagonal(); //V*D^{-1/2}
      //now generate random numbers...
      NumericVector r(iter*pos);
      r = rnorm(iter*pos, 0, 1); // using vectorization from Rcpp sugar
      Map<VectorXd> rvec(as<Map<VectorXd> >(r));
      Map<MatrixXd> rmat(rvec.data(), pos, iter);
      MatrixXd samp(pos, iter);
      // Eigen::LLT<MatrixXd> hesssqrt;
      // hesssqrt.compute(hess);
      // if (hesssqrt.info() != 1){
      //   Rcpp::stop("Cholesky of Hessian failed, probably not positive definite");
      // }
      // NumericVector r(iter*N*(D-1));
      // r = rnorm(iter*N*(D-1), 0, 1); // using vectorization from Rcpp sugar
      // Map<VectorXd> rvec(as<Map<VectorXd> >(r));
      // Map<MatrixXd> rmat(rvec.data(), N*(D-1), iter);
      // MatrixXd samp(N*(D-1), iter);
      // samp = hesssqrt.matrixL().solve(rmat);  // calculate errors of approximation
      samp = hesssqrt*rmat;
      samp.colwise() += eta; // add mean of approximation
      IntegerVector d = IntegerVector::create(D-1, N, iter);
      NumericVector samples = wrap(samp);
      samples.attr("dim") = d; // convert to 3d array for return to R
      out[4] = samples;
    }
  }
  return out;
}