#include <PibbleCollapsed.h>
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;


//' Calculations for the Collapsed Pibble Model
//'
//' Functions providing access to the Log Likelihood, Gradient, and Hessian
//' of the collapsed pibble model. Note: These are convenience functions
//' but are not as optimized as direct coding of the PibbleCollapsed
//' C++ class due to a lack of Memoization. By contrast function optimPibbleCollapsed
//' is much more optimized and massively cuts down on repeated calculations.
//' A more efficient Rcpp module based implementation of these functions
//' may following if the future. For model details see \code{\link{optimPibbleCollapsed}}
//' documentation
//'
//' @inheritParams optimPibbleCollapsed
//' @param AInv Inverse of A for LTP (for Pibble defined as 
//'   AInv = solve(diag(N)+ t(X)\%*\%Gamma\%*\%X) )
//' @param KInv Inverse of K for LTP (for Pibble defined as KInv = solve(Xi))
//' @param eta matrix (D-1)xN of parameter values at which to calculate quantities
//' @param sylv (default:false) if true and if N < D-1 will use sylvester determinant
//'   identity to speed computation
//' @return see below
//'   \itemize{
//'     \item loglikPibbleCollapsed - double
//'     \item gradPibbleCollapsed - vector
//'     \item hessPibbleCollapsed- matrix
//'   }
//' @md
//' @export
//' @examples
//' D <- 10
//' Q <- 2
//' N <- 30
//'
//' # Simulate Data
//' Sigma <- diag(sample(1:8, D-1, replace=TRUE))
//' Sigma[2, 3] <- Sigma[3,2] <- -1
//' Gamma <- diag(sqrt(rnorm(Q)^2))
//' Theta <- matrix(0, D-1, Q)
//' Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1)), nrow=D-1)%*%chol(Gamma)
//' X <- matrix(rnorm(N*(Q-1)), Q-1, N)
//' X <- rbind(1, X)
//' Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1)), nrow=D-1)
//' Pi <- t(driver::alrInv(t(Eta)))
//' Y <- matrix(0, D, N)
//' for (i in 1:N) Y[,i] <- rmultinom(1, sample(5000:10000), prob = Pi[,i])
//'
//' # Priors
//' upsilon <- D+10
//' Xi <- Sigma*(upsilon-D)
//'
//' # Precompute
//' KInv <- solve(Xi)
//' AInv <- solve(diag(N)+ t(X)%*%Gamma%*%X)
//' ThetaX <- Theta%*%X
//'
//'
//' loglikPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, Eta)
//' gradPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, Eta)[1:5]
//' hessPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, Eta)[1:5,1:5]
// [[Rcpp::export]]
double loglikPibbleCollapsed(const Eigen::ArrayXXd Y,
                  const double upsilon,
                  const Eigen::MatrixXd ThetaX,
                  const Eigen::MatrixXd KInv,
                  const Eigen::MatrixXd AInv,
                  Eigen::MatrixXd eta,
                  bool sylv=false){
  // note inverting naming structure here to accord with manuscript
  PibbleCollapsed cm(Y, upsilon, ThetaX, KInv, AInv, sylv); 
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec);
  return cm.calcLogLik(etavec);
}


//' @rdname loglikPibbleCollapsed
//' @export
// [[Rcpp::export]]
Eigen::VectorXd gradPibbleCollapsed(const Eigen::ArrayXXd Y,
                         const double upsilon,
                         const Eigen::MatrixXd ThetaX,
                         const Eigen::MatrixXd KInv,
                         const Eigen::MatrixXd AInv,
                         Eigen::MatrixXd eta,
                         bool sylv=false){
  // note inverting naming structure here to accord with manuscript
  PibbleCollapsed cm(Y, upsilon, ThetaX, KInv, AInv, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec);
  cm.updateWithEtaGH();
  return cm.calcGrad();
}

//' @rdname loglikPibbleCollapsed
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd hessPibbleCollapsed(const Eigen::ArrayXXd Y,
                         const double upsilon,
                         const Eigen::MatrixXd ThetaX,
                         const Eigen::MatrixXd KInv,
                         const Eigen::MatrixXd AInv,
                         Eigen::MatrixXd eta,
                         bool sylv=false){
  // note inverting naming structure here to accord with manuscript
  PibbleCollapsed cm(Y, upsilon, ThetaX, KInv, AInv, sylv);
  Map<VectorXd> etavec(eta.data(), eta.size());
  cm.updateWithEtaLL(etavec);
  cm.updateWithEtaGH();
  return cm.calcHess();
}

// //' Hessian Vector Product using Finite Differences
// //' @rdname hessVectorProd
// //' @export
// // [[Rcpp::export]]
// Eigen::VectorXd hessVectorProd(const Eigen::ArrayXXd Y,
//                          const double upsilon,
//                          const Eigen::MatrixXd ThetaX,
//                          const Eigen::MatrixXd KInv,
//                          const Eigen::MatrixXd AInv,
//                          Eigen::MatrixXd eta,
//                          Eigen::VectorXd v,
//                          double r,
//                          bool sylv=false){
//   // note inverting naming structure here to accord with manuscript
//   PibbleCollapsed cm(Y, upsilon, ThetaX, KInv, AInv, sylv);
//   Map<VectorXd> etavec(eta.data(), eta.size());
//   return cm.calcHessVectorProd(etavec, v, r);
// }

// //' Backtracking line search
// //' @rdname lineSearch
// //' @export
// // [[Rcpp::export]]
// Eigen::VectorXd lineSearch(const Eigen::ArrayXXd Y,
//                          const double upsilon,
//                          const Eigen::MatrixXd ThetaX,
//                          const Eigen::MatrixXd K,
//                          const Eigen::MatrixXd A,
//                          Eigen::MatrixXd eta,
//                          int direction,
//                          double rho,
//                          double c) {
//   // direction: element of eta to step gradient along; (+)
//   // rho: the backtrack step between (0,1) usually 0.5
//   // c: parameter between 0 and 1, usually 0.0001
//   // calculate gradient at current eta
//   Eigen::VectorXd grad = gradPibbleCollapsed(Y, upsilon, ThetaX, K, A, eta);
//   // choose d
//   Eigen::VectorXd d = Eigen::VectorXd::Zero(eta.size());
//   // R indexing to C indexing?
//   d(direction-1) = 1;
//   if(grad(direction-1) < 0) {
//     d(direction-1) = -1;
//   }
//   // how to set initial forward step size?
//   double step = 100;
//   double f0 = loglikPibbleCollapsed(Y, upsilon, ThetaX, K, A, eta);
//   // printf("Original likelihood: %f\n", f0);
//   VectorXd new_eta = eta + step*d;
//   double f1 = loglikPibbleCollapsed(Y, upsilon, ThetaX, K, A, new_eta);
//   // printf("New likelihood: %f\n", f1);
//   // we want an increase in llik, hence the stopping condition
//   while (f1 < f0 + c*step*grad.transpose()*d) {
//     step = step*rho;
//     new_eta = eta + step*d;
//     f1 = loglikPibbleCollapsed(Y, upsilon, ThetaX, K, A, new_eta);
//     // printf("New likelihood: %f\n", f1);
//   }
//   // terminates with an optimal step size
//   return (new_eta);
// }
