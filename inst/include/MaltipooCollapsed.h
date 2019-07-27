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
 *  Where A = (I_N + e^{ell_1}*X*U_1*X' + ... + e^{ell_P}*X*U_P*X' )^{-1},
 *  K^{-1} =Xi is a D-1xD-1 covariance 
 *  
 *  Currently treats ell as a fixed parameter to be estimated
 *  by MAP. 
 *  
 *  matrix, and U_1,...U_P are Q x Q covariance matrix
 */
class MaltipooCollapsed : public Numer::MFuncGrad
{
  private:
    const ArrayXXd Y;
    const double upsilon;
    const MatrixXd Theta;
    const MatrixXd X;
    MatrixXd ThetaX;
    const MatrixXd K; // passed as Xi^{-1}
    const MatrixXd U; // PQ x Q matrix of concatenated deltas
    MatrixXd XTUX;
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
    //Eigen::ColPivHouseholderQR<MatrixXd> Sdec;
    //Eigen::ColPivHouseholderQR<MatrixXd> Ainvdec;
    Eigen::PartialPivLU<MatrixXd> Sdec;
    Eigen::PartialPivLU<MatrixXd> Ainvdec;
    MatrixXd E;  // eta-ThetaX
    ArrayXXd O;  // exp{eta}
    // only needed for gradient and hessian
    MatrixXd rhomat;
    VectorXd rho; 
    MatrixXd C;
    MatrixXd R;
    MatrixXd M; // for maltipoo specifically
    bool sylv;
    
  public:
    MaltipooCollapsed(const ArrayXXd Y_,          // constructor
                        const double upsilon_,
                        const MatrixXd Theta_,
                        const MatrixXd X_,
                        const MatrixXd K_,
                        const MatrixXd U_,
                        bool sylv=false) :
    Y(Y_), upsilon(upsilon_), Theta(Theta_), X(X_), K(K_), U(U_)
    {
      D = Y.rows();           // number of multinomial categories
      N = Y.cols();           // number of samples
      Q = X.rows();
      P = U.rows()/Q;
      ThetaX.noalias() = Theta*X;
      n = Y.colwise().sum();  // total number of counts per sample
      delta = 0.5*(upsilon + N + D - 2.0);
      this->sylv = sylv;
      XTUX = MatrixXd::Zero(P*N, N);
      for (int i=0; i<P; i++){
        XTUX.middleRows(N*i, N).noalias() = X.transpose()*U.middleRows(Q*i, Q)*X;
      }
    }
    ~MaltipooCollapsed(){}                      // destructor
    
    // Update with Eta when it comes in as a vector
    void updateWithEtaLL(const Ref<const VectorXd>& etavec, const Ref<const VectorXd>& ell){
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      E = eta - ThetaX;
      
      Ainv = MatrixXd::Identity(N, N);
      for (int i=0; i<P; i++){
        Ainv += exp(ell(i))*XTUX.middleRows(N*i, N);
      }
      Ainvdec.compute(Ainv);
      A = Ainvdec.inverse(); 
      
      if (sylv & (N < (D-1))){
        S.noalias() = A*E.transpose()*K*E;
        S.diagonal() += VectorXd::Ones(N);
      } else {
        S.noalias() = K*E*A*E.transpose();
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
        C.noalias() = K*E;
        R.noalias() = Sdec.solve(A); // S^{-1}AInv
        MatrixXd SinvK = Sdec.solve(K);
        M.noalias() = A*E.transpose()*SinvK*E*A;
      } else {
        C.noalias() = A*E.transpose();
        R.noalias() = Sdec.solve(K); // S^{-1}K
        M.noalias() = A*E.transpose()*R*E*A;
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
      //ll -= 0.5*(D-1)*Ainvdec.logAbsDeterminant();
      ld = 0.0;
      c = Ainvdec.permutationP().determinant();
      diagLU = Ainvdec.matrixLU().diagonal();
      for (unsigned i = 0; i < diagLU.rows(); ++i) {
        const double& lii = diagLU(i);
        if (lii < 0.0) c *= -1;
        ld += log(std::abs(lii));
      }
      ld += log(c);
      ll -= 0.5*(D-1)*ld;
      return ll;
    }
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    VectorXd calcGrad(const Ref<const VectorXd>& ell){ 
      // For Multinomial
      MatrixXd g = (Y.topRows(D-1)  - (rhomat.array().rowwise()*n.array())).matrix();
      // For MatrixVariate T
      if (sylv & (N < (D-1))){
        g.noalias() += -delta*C*(R+R.transpose());
      } else {
        g.noalias() += -delta*(R + R.transpose())*C.transpose();
      }
      Map<VectorXd> eg(g.data(), g.size()); 
      VectorXd sg(P);
      for (int i=0; i<P; i++){
        sg(i) = delta*(M.array()*XTUX.middleRows(N*i, N).array()).sum();
        sg(i) -= 0.5*(D-1)*(A.array()*XTUX.middleRows(N*i,N).array()).sum();
        sg(i) = exp(ell(i))*sg(i);
      }
      VectorXd grad(N*(D-1)+P);
      grad << eg, sg;
      return grad; // not transposing (leaving as vector)
    }
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    MatrixXd calcHess(const Ref<const VectorXd>& ell){
      bool tmp_sylv = sylv;
      if (sylv & (N < (D-1))){
        MatrixXd eta = E + ThetaX;
        Map<VectorXd> etavec(eta.data(), N*(D-1));
        this->sylv=false;
        updateWithEtaLL(etavec, ell);
        updateWithEtaGH();
      }
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
    
    // should return blocks of size D-1 x D-1 stacked in a N(D-1) x D-1 matrix
    MatrixXd calcPartialHess(){
      // For Multinomial only
      MatrixXd H = ArrayXXd::Zero(N*(D-1), D-1);
      MatrixXd W(D-1, D-1);
      VectorXd rhoseg(D-1);
      for (int j=0; j<N; j++){
        rhoseg = rho.segment(j*(D-1), D-1); // rho calculated in updateWithEtaGH()
        W.noalias() = rhoseg*rhoseg.transpose();
        W.diagonal() -= rhoseg;
        H.block(j*(D-1), 0, D-1, D-1).noalias() += n(j)*W;
      }
      return H;
    }
    
    
    // function for use by ADAMOptimizer wrapper (and for RcppNumeric L-BFGS)
    virtual double f_grad(Numer::Constvec& pars, Numer::Refvec grad){
      const Map<const VectorXd> eta(pars.head(N*(D-1)).data(), N*D-1);
      const Map<const VectorXd> ell(pars.tail(P).data(), P);
      updateWithEtaLL(eta, ell);    // precompute things needed for LogLik
      updateWithEtaGH();       // precompute things needed for gradient and hessian
      grad = -calcGrad(ell);      // negative because wraper minimizes
      return -calcLogLik(eta); // negative because wraper minimizes
    }
};


#endif