#include <RcppEigen.h>
#include <MatDist.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;
using Eigen::Map;
using Eigen::Lower;

//Eta should be array with dim [D-1, N, iter]
//' @export
// [[Rcpp::export]]
List uncollapseMMTC(const Eigen::Map<Eigen::VectorXd> eta, // note this is essentially eta
                    const Eigen::Map<Eigen::MatrixXd> X, 
                    const Eigen::Map<Eigen::MatrixXd> Theta,
                    const Eigen::Map<Eigen::MatrixXd> Gamma, 
                    const Eigen::Map<Eigen::MatrixXd> Xi, 
                    const double upsilon){
  List out(2);
  out.names() = CharacterVector::create("Theta", "Sigma");
  int Q = Gamma.rows();
  int D = Xi.rows()+1;
  int N = X.cols();
  int iter = eta.size()/(N*(D-1)); // assumes result is an integer !!!
  double upsilonN = upsilon + N;
  MatrixXd GammaInv(Gamma.lu().inverse());
  MatrixXd GammaInvN(GammaInv + X*X.transpose());
  MatrixXd GammaN(GammaInvN.lu().inverse());
  MatrixXd LGammaN(GammaN.llt().matrixL());
  //const Map<const MatrixXd> Eta(NULL);
  MatrixXd ThetaGammaInvGammaN(Theta*GammaInv*GammaN);
  MatrixXd XGammaN(X.transpose()*GammaN);
  // // Storage for computation
  MatrixXd ThetaN(D-1, Q);
  MatrixXd XiN(D-1, D-1);
  MatrixXd ThetaDraw(D-1, Q);
  MatrixXd LSigmaDraw(D-1, D-1);
  MatrixXd SigmaDraw(D-1, D-1);
  MatrixXd ETheta(D-1, Q);
  MatrixXd EEta(D-1, N);
  // Storage for output
  MatrixXd ThetaDrawO((D-1)*Q, iter);
  MatrixXd SigmaDrawO((D-1)*(D-1), iter);
  // 
  // Rcout << GammaInv << std::endl;
  // Rcout << "  " << std::endl;
  // Rcout << GammaInvN << std::endl;
  // Rcout << "  " << std::endl;
  // Rcout << LGammaN << std::endl;
  // Rcout << "  " << std::endl;
  // Rcout<< ThetaGammaInvGammaN << std::endl;
  // Rcout << "  " << std::endl;
  //Rcout << XGammaN << std::endl;
  
  

  // iterate over all draws of eta
  for (int i=0; i < iter; i++){
    VectorXd EtaV(eta.segment(i*N*(D-1),N*(D-1)));
    Map<MatrixXd> Eta(EtaV.data(), D-1, N);
    //Rcout << Eta.col(1).transpose() << std::endl;
    ThetaN = Eta*XGammaN+ThetaGammaInvGammaN;
    //Rcout << ThetaN.row(1) << std::endl;
    ETheta = ThetaN-Theta;
    EEta = Eta-ThetaN*X;
    XiN = Xi+ EEta*EEta.transpose() + ETheta*Gamma*ETheta.transpose();
    
    // Draw Random Component
    LSigmaDraw = rInvWishCholesky(upsilonN, XiN).matrix();
    ThetaDraw = rMatNormalCholesky(ThetaN, LSigmaDraw, LGammaN.matrix());
    // Rcout << "   " << std::endl;
    // Rcout << LSigmaDraw << std::endl;
    
    // map output to vectors
    Map<VectorXd> ThetaDrawVec(ThetaDraw.data(), ThetaDraw.size());
    ThetaDrawO.col(i) = ThetaDrawVec;
    SigmaDraw = LSigmaDraw*LSigmaDraw.transpose();
    Map<VectorXd> SigmaDrawVec(SigmaDraw.data(), SigmaDraw.size());
    SigmaDrawO.col(i) = SigmaDrawVec;
  }

  IntegerVector dTheta = IntegerVector::create(D-1, Q, iter);
  IntegerVector dSigma = IntegerVector::create(D-1, D-1, iter);
  NumericVector nvTheta = wrap(ThetaDrawO);
  NumericVector nvSigma = wrap(SigmaDrawO);
  nvTheta.attr("dim") = dTheta;
  nvSigma.attr("dim") = dSigma;
  out[0] = nvTheta;
  out[1] = nvSigma;
  return out;
}


// A few functions for testing MatDist Functions
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd rMatNormalCholesky_test(Eigen::MatrixXd M, 
                                        Eigen::MatrixXd LU, 
                                        Eigen::MatrixXd LV){
  return rMatNormalCholesky(M, LU, LV);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd rInvWishCholesky_test(int v, Eigen::MatrixXd Psi){
  return rInvWishCholesky(v, Psi);
}
  
  
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd rMatUnitNormal_test(int n, int m){
  MatrixXd X(n,m);
  fillUnitNormal(X);
  return X;
}