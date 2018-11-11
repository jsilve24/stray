#include <TruncatedLaplaceApproximation.h>
#include <MongrelCollapsed.h>
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using mongreltrunclapap::eigRes;

// A few functions for testing TruncatedLaplaceApproximation.h
// etavec = c(eta) in R
// [[Rcpp::depends(RSpectra)]]
// [[Rcpp::export]]
List MongrelTruncatedEigen_mongrel_test(Eigen::ArrayXXd Y, double upsilon, 
                                Eigen::MatrixXd ThetaX, 
                                Eigen::MatrixXd K, 
                                Eigen::MatrixXd A, 
                                Eigen::VectorXd etavec, 
                                double r, int nev, int ncv){
  MongrelCollapsed cm(Y, upsilon, ThetaX, K, A);
  eigRes eigs = mongreltrunclapap::MongrelTruncatedEigen(cm, etavec, r, nev, ncv);
  List out(2);
  out[0] = eigs.eigenvalues;
  out[1] = eigs.eigenvectors;
  out.names() = CharacterVector::create("values", "vectors");
  return out;
}
