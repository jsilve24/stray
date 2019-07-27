#include <MaltipooCollapsed.h>
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;


//' @rdname loglikMaltipooCollapsed
//' @export
// [[Rcpp::export]]
double loglikMaltipooCollapsed(const Eigen::ArrayXXd Y,
                  const double upsilon,
                  const Eigen::MatrixXd Theta,
                  const Eigen::MatrixXd X,
                  const Eigen::MatrixXd K,
                  const Eigen::MatrixXd U,
                  Eigen::MatrixXd eta,
                  Eigen::VectorXd ell,
                  bool sylv=false){
  MaltipooCollapsed cm(Y, upsilon, Theta, X, K, U, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec, ell);
  return cm.calcLogLik(etavec);
}

//' @rdname gradMaltipooCollapsed
//' @export
// [[Rcpp::export]]
Eigen::VectorXd gradMaltipooCollapsed(const Eigen::ArrayXXd Y,
                         const double upsilon,
                         const Eigen::MatrixXd Theta,
                         const Eigen::MatrixXd X,
                         const Eigen::MatrixXd K,
                         const Eigen::MatrixXd U,
                         Eigen::MatrixXd eta,
                         Eigen::VectorXd ell,
                         bool sylv=false){
  MaltipooCollapsed cm(Y, upsilon, Theta, X, K, U, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec, ell);
  cm.updateWithEtaGH();
  return cm.calcGrad(ell);
}

//' @rdname hessMaltipooCollapsed
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd hessMaltipooCollapsed(const Eigen::ArrayXXd Y,
                         const double upsilon,
                         const Eigen::MatrixXd Theta,
                         const Eigen::MatrixXd X,
                         const Eigen::MatrixXd K,
                         const Eigen::MatrixXd U,
                         Eigen::MatrixXd eta,
                         Eigen::VectorXd ell,
                         bool sylv=false){
  MaltipooCollapsed cm(Y, upsilon, Theta, X, K, U, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec, ell);
  cm.updateWithEtaGH();
  return cm.calcHess(ell);
}
