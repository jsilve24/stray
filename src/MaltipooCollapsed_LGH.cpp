#include <MaltipooCollapsed.h>
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;

//' Calculations for the Collapsed Maltipoo Model
//'
//' Functions providing access to the Log Likelihood, Gradient, and Hessian
//' of the collapsed maltipoo model. Note: These are convenience functions
//' but are not as optimized as direct coding of the MaltipooCollapsed
//' C++ class due to a lack of Memoization. By contrast function optimMaltipooCollapsed
//' is much more optimized and massively cuts down on repeated calculations.
//' A more efficient Rcpp module based implementation of these functions
//' may following if the future. For model details see \code{\link{optimMaltipooCollapsed}}
//' documentation
//' @inheritParams optimMaltipooCollapsed
//' @param eta matrix (D-1)xN of parameter values at which to calculate quantities
//' @param sylv (default:false) if true and if N < D-1 will use sylvester determinant
//'   identity to speed computation
//' @param ell P-vector of scale factors for each variance component (aka VCScale) 
//' @name loglikMaltipooCollapsed
//' @export
// [[Rcpp::export]]
double loglikMaltipooCollapsed(const Eigen::ArrayXXd Y,
                  const double upsilon,
                  const Eigen::MatrixXd Theta,
                  const Eigen::MatrixXd X,
                  const Eigen::MatrixXd KInv,
                  const Eigen::MatrixXd U,
                  Eigen::MatrixXd eta,
                  Eigen::VectorXd ell,
                  bool sylv=false){
  MaltipooCollapsed cm(Y, upsilon, Theta, X, KInv, U, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec, ell);
  return cm.calcLogLik(etavec);
}

//' @rdname loglikMaltipooCollapsed
//' @export
// [[Rcpp::export]]
Eigen::VectorXd gradMaltipooCollapsed(const Eigen::ArrayXXd Y,
                         const double upsilon,
                         const Eigen::MatrixXd Theta,
                         const Eigen::MatrixXd X,
                         const Eigen::MatrixXd KInv,
                         const Eigen::MatrixXd U,
                         Eigen::MatrixXd eta,
                         Eigen::VectorXd ell,
                         bool sylv=false){
  MaltipooCollapsed cm(Y, upsilon, Theta, X, KInv, U, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec, ell);
  cm.updateWithEtaGH();
  return cm.calcGrad(ell);
}

//' @rdname loglikMaltipooCollapsed
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd hessMaltipooCollapsed(const Eigen::ArrayXXd Y,
                         const double upsilon,
                         const Eigen::MatrixXd Theta,
                         const Eigen::MatrixXd X,
                         const Eigen::MatrixXd KInv,
                         const Eigen::MatrixXd U,
                         Eigen::MatrixXd eta,
                         Eigen::VectorXd ell,
                         bool sylv=false){
  MaltipooCollapsed cm(Y, upsilon, Theta, X, KInv, U, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec, ell);
  cm.updateWithEtaGH();
  return cm.calcHess(ell);
}
