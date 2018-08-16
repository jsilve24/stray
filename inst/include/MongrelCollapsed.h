#ifndef MONGREL_MMTC_H
#define MONGREL_MMTC_H

#include <RcppNumerical.h>
#include <MatrixAlgebra.h>
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
 *    Eta ~ T_{D-1, N}(upsilon, Theta*X, K^{-1}, A^{-1})
 *
 *  Where A = (I_N + X*Gamma*X')^{-1}, K^{-1} =Xi is a D-1xD-1 covariance 
 *  matrix, and Gamma is a Q x Q covariance matrix
 */
class MongrelCollapsed : public Numer::MFuncGrad
{
  private:
    const ArrayXXd Y;
    const double upsilon;
    const MatrixXd ThetaX;
    const MatrixXd K;
    const MatrixXd A;
    // computed quantities 
    int D;
    int N;
    double delta;
    Eigen::ArrayXd m;
    Eigen::RowVectorXd n;
    MatrixXd S;  // I_D-1 + KEAE'
    Eigen::ColPivHouseholderQR<MatrixXd> Sdec;
    MatrixXd E;  // eta-ThetaX
    ArrayXXd O;  // exp{eta}
    // only needed for gradient and hessian
    MatrixXd rhomat;
    VectorXd rho; 
    MatrixXd C;
    MatrixXd R;
    
    // temporary or testing
    int t;
    
    
  public:
    MongrelCollapsed(const ArrayXXd Y_,          // constructor
                        const double upsilon_,
                        const MatrixXd ThetaX_,
                        const MatrixXd K_,
                        const MatrixXd A_) :
    Y(Y_), upsilon(upsilon_), ThetaX(ThetaX_), K(K_), A(A_)
    {
      D = Y.rows();           // number of multinomial categories
      N = Y.cols();           // number of samples
      n = Y.colwise().sum();  // total number of counts per sample
      delta = 0.5*(upsilon + N - D - 2.0);
      
      // temporary or testing
      t = 0;
    }
    ~MongrelCollapsed(){}                      // destructor
    
    // Update with Eta when it comes in as a vector
    void updateWithEtaLL(const Ref<const VectorXd>& etavec){
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      E = eta - ThetaX;
      S.noalias() = K*E*A*E.transpose();
      S.diagonal() += VectorXd::Ones(1, D-1);
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
      C.noalias() = A*E.transpose();
      R.noalias() = Sdec.solve(K); // S^{-1}K
    }
    
    // Must have called updateWithEtaLL first 
    double calcLogLik(const Ref<const VectorXd>& etavec){
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      double ll=0.0;
      // start with multinomial ll
      ll += (Y.topRows(D-1)*eta.array()).sum() - n*m.log().matrix();
      // Now compute collapsed prior ll
      ll -= 0.5*(upsilon+N-D-2)*Sdec.logAbsDeterminant();
      return ll;
    }
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    VectorXd calcGrad(){
      // For Multinomial
      MatrixXd g = (Y - (rhomat.array().rowwise()*n.array())).matrix();
      // For MatrixVariate T
      g.noalias() += -delta*(R + R.transpose())*C.transpose();
      Map<VectorXd> grad(g.data(), g.size()); 
      return grad; // not transposing (leaving as vector)
    }
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    VectorXd calcGrad_wnoise(double sigma, double gamma){
      // For Multinomial
      MatrixXd g = (Y - (rhomat.array().rowwise()*n.array())).matrix();
      // For MatrixVariate T
      g.noalias() += -delta*(R + R.transpose())*C.transpose();
      Map<VectorXd> grad(g.data(), g.size()); 
      
      // add noise
      // NumericVector gnoise(N*(D-1));
      //double ss= sigma/pow(1+t, gamma);
      // gnoise = ss*rnorm(N*(D-1), 0, 1);
      // t++;
      // Rcout << t << " ss: " << ss << std::endl;
      // Rcout << "grad norm: " <<grad.norm() << std::endl;
      // Map<VectorXd> noise(as<Map<VectorXd> >(gnoise));
      // Rcout << "noise norm: " << noise.norm() << std::endl;
      // return grad+noise; // not transposing (leaving as vector)
      
      // add noise as fraction of gradient
      NumericVector gnoise(N*(D-1));
      gnoise = rnorm(N*(D-1), 0, 1);
      //double ss= sigma/pow(1+t, gamma);
      double ss = pow(sigma, gamma*(1+t));
      t++;
      Map<Eigen::ArrayXd> noise(as<Map<Eigen::ArrayXd> >(gnoise));
      return grad+(ss*grad.array().abs()*noise).matrix(); // not transposing (leaving as vector)
    }
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    MatrixXd calcHess(){
      // for MatrixVariate T
      MatrixXd H(N*(D-1), N*(D-1));
      MatrixXd RCT(D-1, N);
      MatrixXd CR(N, D-1);
      MatrixXd L(N*(D-1), N*(D-1));
      RCT.noalias() = R*C.transpose();
      CR.noalias() = C*R;
      L.noalias() = krondense(C*RCT, R.transpose());
      H.noalias() = krondense(A, R+R.transpose());
      H.noalias() -= L+L.transpose();
      L.noalias() = krondense(RCT, RCT.transpose()); // reuse L
      L.noalias() += krondense(CR.transpose(), CR); // reuse L
      H.noalias() -= tveclmult(N, D-1, L);
      H.noalias() = -delta * H;
      // For Multinomial
      MatrixXd W(D-1, D-1);
      VectorXd rhoseg(D-1);
      for (int j=0; j<N; j++){
        rhoseg = rho.segment(j*(D-1), D-1);
        W.noalias() = rhoseg*rhoseg.transpose();
        W.diagonal() -= rhoseg;
        H.block(j*(D-1), j*(D-1), D-1, D-1).noalias()  += n(j)*W;
      }
      return H;
    }
    
    // function for use by ADAMOptimizer wrapper (and for RcppNumeric L-BFGS)
    virtual double f_grad(Numer::Constvec& eta, Numer::Refvec grad){
      updateWithEtaLL(eta);    // precompute things needed for LogLik
      updateWithEtaGH();       // precompute things needed for gradient and hessian
      grad = -calcGrad();      // negative because wraper minimizes
      //grad = -calcGrad_wnoise(.3, 1.001); // optional but not very useful
      return -calcLogLik(eta); // negative because wraper minimizes
    }
};


#endif