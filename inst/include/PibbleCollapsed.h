#ifndef MONGREL_MMTC_H
#define MONGREL_MMTC_H

#include <MatrixAlgebra.h>
#include <MongrelModelClass.h>

#ifdef STRAY_USE_MKL
 #include <mkl.h>
#endif 

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::Ref;

/* Class implementing LogLik, Gradient, and Hessian calculations
 *  for the Multinomial Matrix-T collapsed model. 
 *  
 *  Notation: Let Z_j denote the J-th row of a matrix Z.
 *  
 *  Model:
 *    Y_j ~ Multinomial(Pi_j)
 *    Pi_j = Phi^{-1}(Eta_j)   // Phi^{-1} is ALRInv_D transform
 *    Eta ~ T_{D-1, N}(upsilon, Theta*X, K, A)
 *
 *  Where A = (I_N + X*Gamma*X'), K = Xi is a D-1xD-1 covariance 
 *  matrix, and Gamma is a Q x Q covariance matrix
 */
class PibbleCollapsed : public mongrel::MongrelModel {
  private:
    const ArrayXXd Y;
    const double upsilon;
    const MatrixXd ThetaX;
    const MatrixXd KInv;
    const MatrixXd AInv;
    // computed quantities 
    int D;
    int N;
    double delta;
    Eigen::ArrayXd m;
    Eigen::RowVectorXd n;
    MatrixXd S;  // I_D-1 + KEAE'
    //Eigen::HouseholderQR<MatrixXd> Sdec;
    Eigen::PartialPivLU<MatrixXd> Sdec;
    MatrixXd E;  // eta-ThetaX
    ArrayXXd O;  // exp{eta}
    // only needed for gradient and hessian
    MatrixXd rhomat;
    VectorXd rho; 
    MatrixXd C;
    MatrixXd R;
    
    // testing
    bool sylv;
    
    
  public:
    PibbleCollapsed(const ArrayXXd Y_,          // constructor
                        const double upsilon_,
                        const MatrixXd ThetaX_,
                        const MatrixXd KInv_,
                        const MatrixXd AInv_, 
                        bool sylv=false) :
    Y(Y_), upsilon(upsilon_), ThetaX(ThetaX_), KInv(KInv_), AInv(AInv_)
    {
      D = Y.rows();           // number of multinomial categories
      N = Y.cols();           // number of samples
      n = Y.colwise().sum();  // total number of counts per sample
      delta = 0.5*(upsilon + N + D - 2.0);
      this->sylv = sylv;
    }
    ~PibbleCollapsed(){}                      // destructor
    
    // Update with Eta when it comes in as a vector
    void updateWithEtaLL(const Ref<const VectorXd>& etavec){
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      E = eta - ThetaX;
      if (sylv & (N < (D-1))){
        S.noalias() = AInv*E.transpose()*KInv*E;
        S.diagonal() += VectorXd::Ones(N);
      } else {
        S.noalias() = KInv*E*AInv*E.transpose();
        S.diagonal() += VectorXd::Ones(D-1);  
      }
      Sdec.compute(S);
      O = eta.array().exp();
      m = O.colwise().sum();
      m += Eigen::ArrayXd::Ones(N);
    }
    
    // Must be called after updateWithEtaLL 
    void updateWithEtaGH(){
      rhomat = (O.rowwise()/m.transpose()).matrix();
      Map<VectorXd> rhovec(rhomat.data() , rhomat.size());
      rho = rhovec; // probably could be done in one line rather than 2 (above)
      if (sylv & (N < (D-1))){
        C.noalias() = KInv*E;
        R.noalias() = Sdec.solve(AInv); // S^{-1}AInv
      } else {
        C.noalias() = AInv*E.transpose();
        R.noalias() = Sdec.solve(KInv); // S^{-1}KInv    
      }
    }
    
    // Must have called updateWithEtaLL first 
    double calcLogLik(const Ref<const VectorXd>& etavec){
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      double ll=0.0;
      // start with multinomial ll
      ll += (Y.topRows(D-1)*eta.array()).sum() - n*m.log().matrix();
      // Now compute collapsed prior ll
      //ll -= delta*Sdec.logAbsDeterminant();
      // Following was adapted from : 
      //   https://gist.github.com/redpony/fc8a0db6b20f7b1a3f23
      double ld = 0.0;
      double c = Sdec.permutationP().determinant();
      VectorXd diagLU = Sdec.matrixLU().diagonal();
      for (unsigned i = 0; i < diagLU.rows(); ++i) {
        const double& lii = diagLU(i);
        if (lii < 0.0) c *= -1;
        ld += log(std::abs(lii));
      }
      ld += log(c);
      ll -= delta*ld;
      return ll;
    }
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    VectorXd calcGrad(){
      // For Multinomial
      MatrixXd g = (Y.topRows(D-1) - (rhomat.array().rowwise()*n.array())).matrix();
      //Rcout << "dim Y:" << Y.size() << std::endl;
      //Rcout << "dim g multinomial: " << g.size() << std::endl;
      //Rcout << "dim g t: " << (delta*C*(R+R.transpose()).eval()).size() << std::endl;
      // For MatrixVariate T
      if (sylv & (N < (D-1))){
        g.noalias() += -delta*C*(R+R.transpose());
      } else {
        g.noalias() += -delta*(R + R.transpose())*C.transpose();        
      }
      Map<VectorXd> grad(g.data(), g.size()); 
      return grad; // not transposing (leaving as vector)
    }
    
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    MatrixXd calcHess(){
      bool tmp_sylv = sylv;
      if (sylv & (N < (D-1))){
        MatrixXd eta = E + ThetaX;
        Map<VectorXd> etavec(eta.data(), N*(D-1));
        this->sylv=false;
        updateWithEtaLL(etavec);
        updateWithEtaGH();
      }
      // for MatrixVariate T
      MatrixXd H(N*(D-1), N*(D-1));
      MatrixXd RCT(D-1, N);
      MatrixXd CR(N, D-1);
      MatrixXd L(N*(D-1), N*(D-1));
      RCT.noalias() = R*C.transpose();
      CR.noalias() = C*R;
      krondense_inplace(L, C*RCT, R.transpose());
      //L.noalias() = krondense(C*RCT, R.transpose());
      krondense_inplace(H, AInv, R+R.transpose());
      //H.noalias() = krondense(AInv, R+R.transpose());
      H.noalias() -= L+L.transpose();
      krondense_inplace(L, RCT, RCT.transpose());
      //L.noalias() = krondense(RCT, RCT.transpose()); // reuse L
      krondense_inplace_add(L, CR.transpose(), CR);
      //L.noalias() += krondense(CR.transpose(), CR); // reuse L
      tveclmult_minus(N, D-1, L, H);
      //H.noalias() -= tveclmult(N, D-1, L);
      H.noalias() = -delta * H;
      
      // For Multinomial
      VectorXd rho_parallel;
      VectorXd n_parallel;
      rho_parallel = rho; 
      n_parallel = n;
      
      #pragma omp parallel shared(rho_parallel, n_parallel)
      {
      MatrixXd W(D-1, D-1);
      //VectorXd rhoseg(D-1);
      #pragma omp for 
      for (int j=0; j<N; j++){
        //rhoseg = rho.segment(j*(D-1), D-1);
        Eigen::Ref<VectorXd> rhoseg = rho_parallel.segment(j*(D-1), D-1);
        W.noalias() = rhoseg*rhoseg.transpose();
        W.diagonal() -= rhoseg;
        H.block(j*(D-1), j*(D-1), D-1, D-1).noalias()  += n_parallel(j)*W;
      }
      }
      // Turn back on sylv option if it was wanted:
      this->sylv = tmp_sylv;
      return H;
    }
    
    // function to quickly calculate approximation of hessian-vector product
    //  @param etavec eta at which to calculate hessian
    //  @param v vector to multiply by
    //  @param r size of hessian-vector product difference 
    //  @ref https://justindomke.wordpress.com/2009/01/17/hessian-vector-products/
    VectorXd calcHessVectorProd(const Ref<const VectorXd>& etavec, 
                                VectorXd v, double r=0.001){
      updateWithEtaLL(etavec+r*v);
      updateWithEtaGH();
      VectorXd g1 = calcGrad();
      updateWithEtaLL(etavec-r*v);
      updateWithEtaGH();
      VectorXd g2 = calcGrad();
      return (g1-g2).array()/(2.0*r);
    }
    
    int getN() { return N; }
    int getD() { return D; }
    
    
    // function for use by ADAMOptimizer wrapper (and for RcppNumeric L-BFGS)
    virtual double f_grad(Numer::Constvec& eta, Numer::Refvec grad){
      updateWithEtaLL(eta);    // precompute things needed for LogLik
      updateWithEtaGH();       // precompute things needed for gradient and hessian
      grad = -calcGrad();      // negative because wraper minimizes
      return -calcLogLik(eta); // negative because wraper minimizes
    }
    
};






#endif
