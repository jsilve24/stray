#include <RcppEigen.h>
#include "MatrixAlgebra.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Ref;

// computes L %x% R for dense L and R
MatrixXd krondense(const Ref<const MatrixXd>& L, const Ref<const MatrixXd>& R){
  int lr = L.rows();
  int lc = L.cols();
  int rr = R.rows();
  int rc = R.cols();
  MatrixXd out(lr*rr, lc*rc);
  
  for (int i=0; i < lr; i++){
    for (int j=0; j < lc; j++){
      out.block(i*rr, j*rc, rr, rc) = L(i,j)*R;
    }
  }
  return out;
}

// computes TVEC(m,n)*A for mxn matrix A
MatrixXd tveclmult(const int m, const int n, const Ref<const MatrixXd>& A){
  int ar = A.rows();
  int ac = A.cols();
  MatrixXd out(ar, ac);
  
  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      out.row(i*n+j) = A.row(j*m+i);
    }
  }
  return out;
}

