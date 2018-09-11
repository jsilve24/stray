#ifndef MALTIPOO_MMTC_H
#define MALTIPOO_MMTC_H

#include <RcppNumerical.h>
#include <MatrixAlgebra.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::Ref;

/* Class implementing LogLik, Gradient, and Hessian calculations
 *  for the Multinomial Matrix-T collapsed model with Variance Components. 
 *  
 *  Notation: Let Z_j denote the J-th row of a matrix Z.
 *  
 *  Model:
 *    Y_j ~ Multinomial(Pi_j)
 *    Pi_j = Phi^{-1}(Eta_j)   // Phi^{-1} is ALRInv_D transform
 *    Eta ~ T_{D-1, N}(upsilon, Theta*X, K^{-1}, A^{-1})
 *
 *  Where A = (I_N + sigma^2_1*X*Sigma_1*X' + ... + sigma^2_P*X*Sigma_P*X' )^{-1},
 *  K^{-1} =Xi is a D-1xD-1 covariance 
 *  
 *  Currently treats sigma (little sigma) as a fixed parameter to be estimated
 *  by MAP. 
 *  
 *  matrix, and Sigma_1,...Sigma_P are Q x Q covariance matrix
 */
class MaltipooCollapsed : public Numer::MFuncGrad
{
  private:
    const ArrayXXd Y;
    const double upsilon;
    const MatrixXd Theta;
    const MatrixXd X;
    MatrixXd ThetaX;
    const MatrixXd K;
    const MatrixXd Sigma; // PQ x Q matrix of concatenated sigmas
    MatrixXd XTSigmaX;
    MatrixXd A; // no longer constant
    MatrixXd Ainv; // no longer constant
    // computed quantities 
    int D;
    int N;
    int P;
    int Q;
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
    MatrixXd M; // for maltipoo specifically
    
    
  public:
    MaltipooCollapsed(const ArrayXXd Y_,          // constructor
                        const double upsilon_,
                        const MatrixXd Theta_,
                        const MatrixXd X_,
                        const MatrixXd K_,
                        const MatrixXd Sigma_) :
    Y(Y_), upsilon(upsilon_), Theta(Theta_), X(X_), K(K_), Sigma(Sigma_)
    {
      D = Y.rows();           // number of multinomial categories
      N = Y.cols();           // number of samples
      Q = X.rows();
      P = Sigma.rows()/Q;
      ThetaX.noalias() = Theta*X;
      n = Y.colwise().sum();  // total number of counts per sample
      delta = 0.5*(upsilon + N - D - 2.0);
      XTSigmaX = MatrixXd::Zero(P*N, N);
      for (int i=0; i<P; i++){
        XTSigmaX.middleRows(N*i, N).noalias() = X.transpose()*Sigma.middleRows(Q*i, Q)*X;
      }
    }
    ~MaltipooCollapsed(){}                      // destructor
    
    // Update with Eta when it comes in as a vector
    void updateWithEtaLL(const Ref<const VectorXd>& etavec, const Ref<const VectorXd>& sigmavec){
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      E = eta - ThetaX;
      
      Ainv = MatrixXd::Identity(N, N);
      for (int i=0; i<P; i++){
        Ainv += sigmavec(i)*XTSigmaX.middleRows(N*i, N);
      }
      //Eigen::FullPivLU<MatrixXd> lu(Ainv);
      A = Ainv.lu().inverse(); 
        
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
      M.noalias() = Ainv*E.transpose()*R*E*Ainv;
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
      Map<VectorXd> eg(g.data(), g.size()); 
      VectorXd sg(P);
      for (int i=0; i<P; i++){
        sg(i)=-(M.array()*XTSigmaX.middleRows(N*i, N).array()).sum();
      }
      VectorXd grad(N*(D-1)+P);
      grad << eg, sg;
      return grad; // not transposing (leaving as vector)
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
    virtual double f_grad(Numer::Constvec& pars, Numer::Refvec grad){
      const Map<const VectorXd> eta(pars.head(N*(D-1)).data(), N*D-1);
      const Map<const VectorXd> sigma(pars.tail(P).data(), P);
      updateWithEtaLL(eta, sigma);    // precompute things needed for LogLik
      updateWithEtaGH();       // precompute things needed for gradient and hessian
      grad = -calcGrad();      // negative because wraper minimizes
      return -calcLogLik(eta); // negative because wraper minimizes
    }
};


#endif