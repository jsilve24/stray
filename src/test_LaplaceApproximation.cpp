#include <LaplaceApproximation.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// A few functions for testing LaplaceApproximation.h
// [[Rcpp::export]]
Eigen::MatrixXd eigen_lap_test(int n_samples, Eigen::VectorXd m, 
                               Eigen::MatrixXd S, double eigvalthresh){
  lapap::lappars pars = lapap::init_lappars(eigvalthresh);
  int p = m.rows();
  MatrixXd z = MatrixXd::Zero(p, n_samples);
  int status = lapap::eigen_lap(z, m, S, pars);
  if (status==1) Rcpp::stop("decomposition failed");
  return z;
}

// [[Rcpp::export]]
Eigen::MatrixXd cholesky_lap_test(int n_samples, Eigen::VectorXd m, 
                               Eigen::MatrixXd S, double eigvalthresh){
  lapap::lappars pars = lapap::init_lappars(eigvalthresh);
  int p = m.rows();
  MatrixXd z = MatrixXd::Zero(p, n_samples);
  int status = lapap::cholesky_lap(z, m, S, pars);
  if (status==1) Rcpp::stop("decomposition failed");
  return z;
}



// [[Rcpp::export]]
Eigen::MatrixXd LaplaceApproximation_test(int n_samples, Eigen::VectorXd m, 
                                          Eigen::MatrixXd S, String decomp_method, 
                                          double eigvalthresh){
  int p=m.rows();
  MatrixXd z = MatrixXd::Zero(p, n_samples);
  double logInvNegHessDet;
  int status = lapap::LaplaceApproximation(z, m, S, decomp_method, 
                                           eigvalthresh,0, 
                                           logInvNegHessDet);
  if (status==1) Rcpp::stop("decomposition failed");
  return z;
}