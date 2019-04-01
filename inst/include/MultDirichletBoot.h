#ifndef MONGREL_MULTDIRICHLETBOOT_H
#define MONGREL_MULTDIRICHLETBOOT_H

#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::Ref;

namespace MultDirichletBoot{
  
  template <typename T>
  MatrixXd alrInv_default(Eigen::MatrixBase<T>& eta){
    int D = eta.rows()+1;
    int N = eta.cols();
    MatrixXd pi = MatrixXd::Zero(D, N);
    pi.topRows(D-1) = eta;
    pi.array() = pi.array().exp();
    pi.array().rowwise() /= pi.colwise().sum().array();
    return pi;
  }
  
  template <typename T>
  MatrixXd alr_default(Eigen::MatrixBase<T>& pi){
    int D = pi.rows();
    int N = pi.cols();
    MatrixXd eta(D-1,N);
    eta = pi.topRows(D-1);
    eta.array().rowwise() /= pi.row(D-1).array();
    return eta.array().log();
  }
  
  // Sample dirichlet - alpha must be a vector
  template <typename Derived>
  MatrixXd rDirichlet(int n_samples, Eigen::MatrixBase<Derived>& alpha){
    int D = alpha.rows();
    int p = alpha.cols();
    if (p > 1) Rcpp::stop("rDirichlet must only be passed alpha as a vector");
    NumericVector r(n_samples);
    MatrixXd s(D, n_samples);
    for (int i=0; i<D; i++){
      r = rgamma(n_samples, alpha(i), 1);
      Map<VectorXd> rvec(as<Map<VectorXd> >(r));
      s.row(i) = rvec.transpose();
    }
    s.array().rowwise() /= s.colwise().sum().array();
    return s;
  }
  
  
  template <typename T1>
  MatrixXd MultDirichletBoot(int n_samples, Eigen::MatrixBase<T1>& eta, 
                             ArrayXXd Y, double pseudocount){
    int D = eta.rows()+1;
    int N = eta.cols();
    MatrixXd alpha = alrInv_default(eta);
    alpha.array().rowwise() *= Y.colwise().sum();
    alpha.array() += pseudocount; 
    MatrixXd samp(N*(D-1), n_samples);
    MatrixXd s(D, n_samples);
    VectorXd a;
    for (int i=0; i<N; i++){
      a = alpha.col(i);
      s = rDirichlet(n_samples, a);
      // transform to eta
      samp.middleRows(i*(D-1), D-1) = alr_default(s);
    }
    return samp;
  }
}


#endif