#include <RcppEigen.h>
#include <MatDist.h>
using namespace Rcpp;

// [[Rcpp::export]]
void fillUnitNormal_test(Eigen::Map<Eigen::MatrixXd>& Z){
  fillUnitNormal(Z);
}