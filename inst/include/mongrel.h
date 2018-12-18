// File used to import others in order

#ifdef MONGREL_USE_MKL
  #define EIGEN_USE_MKL_ALL
#else
  #ifdef _OPENMP
    #include <omp.h>
  #endif
#endif 
// [[Rcpp::depends(RcppEigen)]]


#include "MatrixAlgebra.h"
#include "MatDist.h"
#include "MultDirichletBoot.h"
#include "SpecialFunctions.h"
#include "LaplaceApproximation.h"
#include "MongrelCollapsed.h"
#include "MaltipooCollapsed.h"
#include "AdamOptim.h"
#include "TruncatedLaplaceApproximation.h"