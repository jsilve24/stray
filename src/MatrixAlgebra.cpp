#include <RcppEigen.h>
#include "MatrixAlgebra.h"
#ifdef STRAY_USE_MKL
  #include <mkl.h>
#endif

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
  
  #pragma omp parallel for
  for (int i=0; i < lr; i++){
    for (int j=0; j < lc; j++){
      out.block(i*rr, j*rc, rr, rc) = L(i,j)*R;
    }
  }
  return out;
}

// computes A=L%x%R (overwrites A)
void krondense_inplace(Ref<MatrixXd> A, const Ref<const MatrixXd>& L, 
                       const Ref<const MatrixXd>& R){
  int lr = L.rows();
  int lc = L.cols();
  int rr = R.rows();
  int rc = R.cols();

  #pragma omp parallel for shared(R, A)
    for (int i=0; i < lr; i++){
      for (int j=0; j < lc; j++){
        A.block(i*rr, j*rc, rr, rc) = L(i,j)*R;
      }
    }
}

// computes A+=L%x%R (overwrites A)
void krondense_inplace_add(Ref<MatrixXd> A, const Ref<const MatrixXd>& L, 
                       const Ref<const MatrixXd>& R){
  int lr = L.rows();
  int lc = L.cols();
  int rr = R.rows();
  int rc = R.cols();
  
#pragma omp parallel for shared(R, A)
  for (int i=0; i < lr; i++){
    for (int j=0; j < lc; j++){
      A.block(i*rr, j*rc, rr, rc) += L(i,j)*R;
    }
  }
}


// computes TVEC(m,n)*A for mxn matrix A
MatrixXd tveclmult(const int m, const int n, const Ref<const MatrixXd>& A){
  int ar = A.rows();
  int ac = A.cols();
  MatrixXd out(ar, ac);
  
  #pragma omp parallel for shared(out, A)
  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++){
      out.row(i*n+j) = A.row(j*m+i);
    }
  }
  return out;
}

// Computes B-TVEC(m,n)*A and overwrites B with result Note:also modifies A!
// assumes A is in column major layout!
void tveclmult_minus(const int m, const int n, Ref<MatrixXd> A, 
                     Ref<MatrixXd> B){
  int ar=A.rows();
  int ac=A.cols();
  
  #if defined(STRAY_USE_MKL)
  Eigen::VectorXi k(ar);
    for (int i=0; i<m; i++){
      for (int j=0; j<n; j++)
        k(i*n+j) = j*m+i+1;
    }
    LAPACKE_dlapmr(LAPACK_COL_MAJOR, true, ar, ac, A.data(), ar, k.data());
    B-=A;
  #else
    #pragma omp parallel for shared(B, A)
    for (int i=0; i<m; i++){
      for (int j=0; j<n; j++)
        B.row(i*n+j) -= A.row(j*m+i);
    }
  #endif 
}

