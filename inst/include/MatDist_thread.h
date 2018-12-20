#ifndef MONGREL_MATDIST_THREAD_H
#define MONGREL_MATDIST_THREAD_H

#include <RcppEigen.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/chi_squared_distribution.hpp>

using Eigen::MatrixXd;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::Lower;


// fills passed dense objects with unit normal random variables
template <typename Derived, typename RNG>
inline void fillUnitNormal_thread(Eigen::DenseBase<Derived>& Z, RNG& rng){
  //trng::normal_dist<> rnorm(0,1);
  boost::random::normal_distribution<> rnorm(0,1);
  int m = Z.rows();
  int n = Z.cols();
  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++)
      Z(i,j) = rnorm(rng);
  }
}


//' The Matrix Normal Distribution (one at a time)
//'
//' Density, and random generation from matrix normal distribution.
//' X ~ MVN(M, U, V) where X is a NxP matrix, U is a NxN Positive Semi-definite
//' covariance matrix over the rows of X and V is a PxP Positive Semi-definite
//' covariance matrix over the columns of X.
//' For Internal use by other C++ Functions
//'
//' @param X NxP Matrix at which to evaluate density
//' @param M NxP Mean Matrix
//' @param LU Lower Cholesky Factor of U such that L_U*t(L_U) = U
//' @param LV Lower Cholesky Factor of V such that L_V*t(L_V) = V
//'
//' @return Matrix
//'
//' @name MatNormal
//' @examples
//' M <- matrix(0, nrow=2, ncol=2)
//' U <- diag(2)
//' V <- matrix(c(1,.5,.5,1), ncol=2)
//' LV <- t(chol(V))
//' LU <- t(chol(U))
//' rMatNormalCholesky(M, LU, LV)
template <typename RNG>
inline Eigen::MatrixXd rMatNormalCholesky_thread(const Eigen::Ref<const Eigen::MatrixXd>& M,
                                          const Eigen::Ref<const Eigen::MatrixXd>& LU,
                                          const Eigen::Ref<const Eigen::MatrixXd>& LV,
                                          RNG& rng){
  int nrows = M.rows();
  int ncols = M.cols();
  MatrixXd Z(nrows, ncols);
  MatrixXd X(nrows, ncols);

  fillUnitNormal_thread(Z, rng);
  X.noalias() =  M + LU*Z*LV.transpose();
  return X;
}

template <typename T, typename RNG>
inline void rMatNormalCholesky_thread_inplace(Eigen::MatrixBase<T>& A, 
                                                 const Eigen::Ref<const Eigen::MatrixXd>& M,
                                                 const Eigen::Ref<const Eigen::MatrixXd>& LU,
                                                 const Eigen::Ref<const Eigen::MatrixXd>& LV,
                                                 RNG& rng){
  int nrows = M.rows();
  int ncols = M.cols();
  MatrixXd Z(nrows, ncols);
  fillUnitNormal_thread(Z, rng);
  A.noalias() =  M + LU*Z*LV.transpose();
}


//' Inverse Wishart Distribution (Wikipedia Parameterization one at a time)
//'
//' Density, and random generation from Inverse Wishart Distribution.
//' W^{-1}(Psi, v) where Psi is a PxP positive definite scale matrix
//' and v is the degrees of freedom with requirement that v > P-1.
//'
//' @param X is positive definite PxP covariance matrix at
//' which to evaluate density
//' @param v degres of freedom (req: v > P-1)
//' @param Psi PxP positive definite scale matrix (covariance matrix)
//'
//' @details Generate Draws from an Inverse Wishart Distribution
//' via the Bartlett Decomposition. Mean is given by Psi/(v-P-1).
//' Mode is given by Psi/(v+P+1). Does minor imput validation to ensure v > P-1,
//' throws range_error if not false.
//'
//' @name InvWishart
//'
//' @return Reverse cholesky factor of the inverse wishart sample. That is it
//' returns an upper triangular matrix U such if V=UU^T, V ~ IW(v, Psi).
template <typename RNG>
inline Eigen::MatrixXd rInvWishRevCholesky_thread(const int v,
                                           const Eigen::Ref<const Eigen::MatrixXd>& Psi,
                                           RNG& rng){
  int p = Psi.rows();
  MatrixXd PsiInv = Psi.llt().solve(MatrixXd::Identity(p,p));
  if (v <= p-1)
    Rcpp::stop("v must be > Psi.rows - 1");
  VectorXd z(p*(p-1)/2);
  fillUnitNormal_thread(z, rng);
  MatrixXd X = MatrixXd::Zero(p, p);
  for (int i=0; i<p; i++){
    //trng::chi_square_dist<> rchisq(v-i);
    boost::random::chi_squared_distribution<> rchisq(v-i);
    X(i,i) = sqrt(rchisq(rng));
    //X(i,i) = sqrt(R::rchisq(v-i)); // zero indexing
  }
  int pos = 0;
  for (int i=1; i<p; i++){
    for (int j=0; j<i; j++){
      X(i,j) = z(pos);
      pos++;
    }
  }
  MatrixXd Y;
  Y.noalias() = PsiInv.llt().matrixL()*X;
  return Y.triangularView<Lower>().solve(MatrixXd::Identity(p,p)).transpose();
}

template <typename T, typename RNG>
inline void rInvWishRevCholesky_thread_inplace(Eigen::MatrixBase<T>& A, 
                                                  const int v,
                                                  const Eigen::Ref<const Eigen::MatrixXd>& Psi,
                                                  RNG& rng){
  int p = Psi.rows();
  MatrixXd PsiInv = Psi.llt().solve(MatrixXd::Identity(p,p));
  if (v <= p-1)
    Rcpp::stop("v must be > Psi.rows - 1");
  VectorXd z(p*(p-1)/2);
  fillUnitNormal_thread(z, rng);
  MatrixXd X = MatrixXd::Zero(p, p);
  for (int i=0; i<p; i++){
    boost::random::chi_squared_distribution<> rchisq(v-i);
    X(i,i) = sqrt(rchisq(rng));
  }
  int pos = 0;
  for (int i=1; i<p; i++){
    for (int j=0; j<i; j++){
      X(i,j) = z(pos);
      pos++;
    }
  }
  A.template noalias() = PsiInv.llt().matrixL()*X;
  A.template triangularView<Lower>().solveInPlace(MatrixXd::Identity(p,p)); // This line is causing problems
  A.template transposeInPlace();
}



#endif
