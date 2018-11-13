#ifndef MONGREL_TRUNCLA_H
#define MONGREL_TRUNCLA_H

#include <MongrelModelClass.h>
#include <SymEigs.h> // from spectra library

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Ref;
using Eigen::Map;
using mongrel::MongrelModel;
using Spectra::SymEigsSolver;
using Spectra::SMALLEST_MAGN;
using Spectra::LARGEST_MAGN;
using Spectra::SMALLEST_ALGE;

namespace mongreltrunclapap{

  struct eigRes {
    int k;
    VectorXd eigenvalues;
    MatrixXd eigenvectors;
  };
  

  // class for Hessian Vector Product matching Spectra Interface
  // NOTE: ADDS NEGATIVE SIGN SO USING NEGATIVE HESSIAN!
  template <typename Scalar>
  class MongrelHVP {
    private:
      MongrelModel& cm;
      const VectorXd etavec;
      double r;
      bool useOffset;
      double offset;
      
  
    public:
      MongrelHVP(MongrelModel& cm_, const VectorXd etavec_, double r_) :  
        cm(cm_), etavec(etavec_), r(r_) {
        this->offset=0.0;
        this->useOffset=false;
      }
      ~MongrelHVP(){}
      
      int rows() const { return cm.getN()*(cm.getD()-1); }
      int cols() const { return cm.getN()*(cm.getD()-1); }
      
      void addOffset(double offset){
        this->offset = offset;
        this->useOffset=true;
      };
      
      void perform_op(const Scalar* x_in, Scalar* y_out) const {
          Map<const VectorXd> x(x_in, cols());
          Map<VectorXd> y(y_out, rows());
          if (useOffset){
            y.noalias() = cm.calcHessVectorProd(etavec, x, r);  
            y.array() -= offset*x.array(); // (H-offset*I) * x
          } else {
            y.noalias() = cm.calcHessVectorProd(etavec, x, r);
          }
      }
  };  
  
  // Simple Definition to Ease Notation Below
  typedef SymEigsSolver<double, LARGEST_MAGN, MongrelHVP<double> > MongrelSymEig;
  
  // Return eigRes object (NOTE: Calculated from negative hessian!)
  // @param cm object inheriting from class mongrel::MongrelModel
  // @param etavec constant vector (MAP estimate really)
  // @param r double used to define finite differences approximation
  // @param nev number of eigenvalues to return (selects smallest magnitude)
  // @param Parameter that controls the convergence speed of the algorithm. 
  //   Typically a larger ncv means faster convergence, but it may also result 
  //   in greater memory use and more matrix operations in each iteration. 
  //   This parameter must satisfy nev<ncv≤n, and is advised to take ncv≥2⋅nev.
  //   (taken directly from SymEigsSolver documentation)
  inline eigRes MongrelTruncatedEigen(MongrelModel& cm, const VectorXd etavec, 
                                        double r, int nev, int ncv){
    // Get Largest Eigenvalue
    MongrelHVP<double> op(cm, etavec, r);
    MongrelSymEig eig_largest(&op, 1, 5);
    eig_largest.init();
    int nconv1 = eig_largest.compute();
    if (eig_largest.info() == Spectra::NOT_CONVERGING){
      Rcpp::stop("Spectra not converging in truncated eigendecomposition");
    }
    if (eig_largest.info() != Spectra::SUCCESSFUL) {
      Rcpp::stop("Something unexpected went wrong in truncated eigendecomposition");
    }
    // 
    // Rcpp::Rcout << "Took " << nconv1 << " iterations to converge First TruncSVD" 
    //             << std::endl;
    double largest = eig_largest.eigenvalues().value(); 
    //Rcpp::Rcout << "Largest: " << largest << std::endl;
    
    op.addOffset(largest); // add negative 

    // Get Smallest part of spectrum
    MongrelSymEig eigs(&op, nev, ncv);

    eigs.init();
    int nconv2 = eigs.compute();
    // 
    // Rcpp::Rcout << "Took " << nconv2 << " iterations to converge First TruncSVD" 
    //             << std::endl;
    
    if (eigs.info() == Spectra::NOT_CONVERGING){
      Rcpp::stop("Spectra not converging in truncated eigendecomposition");
    }
    if (eigs.info() != Spectra::SUCCESSFUL) {
      Rcpp::stop("Something unexpected went wrong in truncated eigendecomposition");
    }
    eigRes res;
    res.k = nev;
    res.eigenvalues=eigs.eigenvalues();
    res.eigenvalues.array() += largest;
    res.eigenvalues.array() *= -1;
    res.eigenvectors = eigs.eigenvectors();
    return res;
  }
  
  template <typename T1, typename T2>  
  // @param cm object inheriting from class mongrel::MongrelModel
  // @param z an object derived from class MatrixBase to overwrite with samples
  // @param m MAP estimate (as a vector)
  // @param k integer (k eigenvectors with smallest overall magnitude)
  // @param r see MongrelTruncatedEigen
  // @param nev see MongrelTruncatedEigen
  // @param ncv see MongrelTruncatedEigen
  // @param eivalthresh threshold for negative 
  //    eigenvalues dictates clipping vs. stopping behavior
  inline int TruncatedLaplaceApproximation(MongrelModel& cm, 
                                           Eigen::MatrixBase<T1>& z, 
                                           Eigen::MatrixBase<T2>& m, 
                                           int k,
                                           double r, 
                                           int nev, int ncv, 
                                           double eigvalthresh){
    int nc = z.cols();
    if (k > cm.getN()*(cm.getD()-1)){
      Rcpp::stop("k must be smaller than N*(D-1)");
    }
    if (k <=0){
      Rcpp::stop("k must be larger than 0");
    }
    
    // NOTE: automatically uses negative hesssian see MongrelHVP code!
    MongrelSymEig eigs = MongrelTruncatedEigen(cm, m, r, nev, ncv);
    
    VectorXd evalinv(eigs.eigenvalues().array().inverse().matrix());
    int excess=0;
    for (int i=0; i<k; i++){
      if (evalinv(i) < eigvalthresh) {
        excess++;
      }
    }
    if (excess > 0){
      Rcpp::warning("Some eigenvalues are below eigvalthresh");
      Rcpp::Rcout << "Eigenvalues" << evalinv.transpose() << std::endl;
      return 1;
    }
    int pos = 0;
    for (int i=0; i<k; i++){
      if (evalinv(pos) > 0)
        pos++;
    }
    if (pos < k) {
      Rcpp::warning("Some small negative eigenvalues are being chopped");
      Rcpp::Rcout << k-pos << " out of " << k <<
        " passed eigenvalue threshold" << std::endl;
    }
    
    MatrixXd invhesssqrt(k, pos);
    invhesssqrt = eigs.eigenvectors().rightCols(pos)*
      evalinv.tail(pos).cwiseSqrt().asDiagonal(); //V*D^{-1/2}
    typename T1::PlainObject samp(pos, nc);
    fillUnitNormal(samp);
    z.template noalias() = invhesssqrt*samp;
    z.template colwise() += m;
    return 0;
  } // end function
  
  
} // end namespace

#endif 
