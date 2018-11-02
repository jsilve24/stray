#include <MultDirichletBoot.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Wrapper functions 
// [[Rcpp::export]]
Eigen::MatrixXd alrInv_default_test(Eigen::MatrixXd eta){
  return MultDirichletBoot::alrInv_default(eta);
}

// Wrapper functions 
// [[Rcpp::export]]
Eigen::MatrixXd alr_default_test(Eigen::MatrixXd pi){
  return MultDirichletBoot::alr_default(pi);
}

// Wrapper functions 
// [[Rcpp::export]]
Eigen::MatrixXd rDirichlet_test(int n_samples, Eigen::VectorXd alpha){
  return MultDirichletBoot::rDirichlet(n_samples, alpha);
}


// Wrapper functions 
// [[Rcpp::export]]
Eigen::MatrixXd MultDirichletBoot_test(int n_samples, Eigen::MatrixXd eta, 
                                       Eigen::ArrayXXd Y, double pseudocount){
  return MultDirichletBoot::MultDirichletBoot(n_samples, eta, Y, pseudocount);
}