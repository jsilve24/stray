#ifndef MONGREL_MATALG_H
#define MONGREL_MATALG_H

#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Ref;

MatrixXd krondense(const Ref<const MatrixXd>& L, const Ref<const MatrixXd>& R);
void krondense_inplace(Ref<MatrixXd> A, const Ref<const MatrixXd>& L, 
                       const Ref<const MatrixXd>& R);
void krondense_inplace_add(Ref<MatrixXd> A, const Ref<const MatrixXd>& L, 
                       const Ref<const MatrixXd>& R);
MatrixXd tveclmult(const int m, const int n, const Ref<const MatrixXd>& A);
void tveclmult_minus(const int m, const int n, Ref<MatrixXd> A, Ref<MatrixXd> B);

#endif