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


//' Uncollapse output from optimMongrelCollapsed to full Mongrel Model
//' 
//' See details for model. Should likely be called following 
//' \code{\link{optimMongrelCollapsed}}. Notation: \code{N} is number of samples,
//' \code{D} is number of multinomial categories, \code{Q} is number
//' of covariates, \code{iter} is the number of samples of \code{eta} (e.g., 
//' the parameter \code{n_samples} in the function \code{optimMongrelCollapsed})
//' 
//' @param eta array of dimension (D-1) x N x iter (e.g., \code{Pars} output of 
//'   function optimMongrelCollapsed)
//' @param X matrix of covariates of dimension Q x N
//' @param Theta matrix of prior mean of dimension (D-1) x Q
//' @param Gamma covariance matrix of dimension Q x Q
//' @param Xi covariance matrix of dimension (D-1) x (D-1)
//' @param upsilon scalar (must be > D) degrees of freedom for InvWishart prior
//' @param ret_mean if true then uses posterior mean of Lambda and Sigma 
//'   corresponding to each sample of eta rather than sampling from 
//'   posterior of Lambda and Sigma (useful if Laplace approximation
//'   is not used (or fails) in optimMongrelCollapsed)
//' 
//' @details Notation: Let Z_j denote the J-th row of a matrix Z.
//' While the collapsed model is given by:
//'    \deqn{Y_j ~ Multinomial(Pi_j)}
//'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
//'    \deqn{Eta ~ T_{D-1, N}(upsilon, Theta*X, K^{-1}, A^{-1})}
//' Where A = (I_N + X * Gamma * X')^{-1}, K^{-1} = Xi is a (D-1)x(D-1) covariance 
//' matrix, Gamma is a Q x Q covariance matrix, and Phi^{-1} is ALRInv_D 
//' transform. 
//' 
//' The uncollapsed model (Full Mongrel model) is given by:
//'    \deqn{Y_j ~ Multinomial(Pi_j)}
//'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
//'    \deqn{Eta ~ MN_{D-1 x N}(Lambda*X, Sigma, I_N)}
//'    \deqn{Lambda ~ MN_{D-1 x Q}(Theta, Sigma, Gamma)}
//'    \deqn{Sigma ~ InvWish(upsilon, Xi)}
//' This function provides a means of sampling from the posterior distribution of 
//' \code{Lambda} and \code{Sigma} given posterior samples of \code{Eta} from 
//' the collapsed model. 
//' @return List with components 
//' 1. Lambda Array of dimension (D-1) x Q x iter (posterior samples)
//' 2. Sigma Array of dimension (D-1) x (D-1) x iter (posterior samples)
//' @export
//' @md
//' @seealso \code{\link{optimMongrelCollapsed}}
//' @examples
//' sim <- mongrel_sim()
//' 
//' # Fit model for eta
//' fit <- optimMongrelCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$K, 
//'                              sim$A, random_mongrel_init(sim$Y))  
//' 
//' # Finally obtain samples from Lambda and Sigma
//' fit2 <- uncollapseMongrelCollapsed(fit$Samples, sim$X, sim$Theta, 
//'                                    sim$Gamma, sim$Xi, sim$upsilon)
// [[Rcpp::export]]
List uncollapseMongrelCollapsed(const Eigen::Map<Eigen::VectorXd> eta, // note this is essentially eta
                    const Eigen::Map<Eigen::MatrixXd> X, 
                    const Eigen::Map<Eigen::MatrixXd> Theta,
                    const Eigen::Map<Eigen::MatrixXd> Gamma, 
                    const Eigen::Map<Eigen::MatrixXd> Xi, 
                    const double upsilon, 
                    bool ret_mean = false){
  List out(2);
  out.names() = CharacterVector::create("Lambda", "Sigma");
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
  MatrixXd XTGammaN(X.transpose()*GammaN);
  // // Storage for computation
  MatrixXd LambdaN(D-1, Q);
  MatrixXd XiN(D-1, D-1);
  MatrixXd LambdaDraw(D-1, Q);
  MatrixXd LSigmaDraw(D-1, D-1);
  MatrixXd SigmaDraw(D-1, D-1);
  MatrixXd ELambda(D-1, Q);
  MatrixXd EEta(D-1, N);
  // Storage for output
  MatrixXd LambdaDrawO((D-1)*Q, iter);
  MatrixXd SigmaDrawO((D-1)*(D-1), iter);

  // iterate over all draws of eta
  for (int i=0; i < iter; i++){
    R_CheckUserInterrupt();
    VectorXd EtaV(eta.segment(i*N*(D-1),N*(D-1)));
    Map<MatrixXd> Eta(EtaV.data(), D-1, N);
    //Rcout << Eta.col(1).transpose() << std::endl;
    LambdaN = Eta*XTGammaN+ThetaGammaInvGammaN;
    //Rcout << LambdaN.row(1) << std::endl;
    ELambda = LambdaN-Theta;
    EEta = Eta-LambdaN*X;
    XiN = Xi+ EEta*EEta.transpose() + ELambda*GammaInv*ELambda.transpose();
    
    if (ret_mean){
      Map<VectorXd> LambdaNVec(LambdaN.data(), LambdaN.size());
      Map<VectorXd> XiNVec(XiN.data(), XiN.size());
      LambdaDrawO.col(i) = LambdaNVec;
      SigmaDrawO.col(i) = (upsilonN-D)*XiNVec; // mean of inverse wishart
    } else {
      // Draw Random Component
      LSigmaDraw = rInvWishRevCholesky(upsilonN, XiN).matrix();
      // Note: Below is valid even though LSigmaDraw is reverse cholesky factor
      LambdaDraw = rMatNormalCholesky(LambdaN, LSigmaDraw, LGammaN.matrix());
      
      // map output to vectors
      Map<VectorXd> LambdaDrawVec(LambdaDraw.data(), LambdaDraw.size());
      LambdaDrawO.col(i) = LambdaDrawVec;
      SigmaDraw = LSigmaDraw*LSigmaDraw.transpose();
      Map<VectorXd> SigmaDrawVec(SigmaDraw.data(), SigmaDraw.size());
      SigmaDrawO.col(i) = SigmaDrawVec;
    }
  }

  IntegerVector dLambda = IntegerVector::create(D-1, Q, iter);
  IntegerVector dSigma = IntegerVector::create(D-1, D-1, iter);
  NumericVector nvLambda = wrap(LambdaDrawO);
  NumericVector nvSigma = wrap(SigmaDrawO);
  nvLambda.attr("dim") = dLambda;
  nvSigma.attr("dim") = dSigma;
  out[0] = nvLambda;
  out[1] = nvSigma;
  
  return out;
}

// A few functions for testing MatDist Functions
// [[Rcpp::export]]
Eigen::MatrixXd rMatNormalCholesky_test(Eigen::MatrixXd M, 
                                        Eigen::MatrixXd LU, 
                                        Eigen::MatrixXd LV){
  return rMatNormalCholesky(M, LU, LV);
}

// [[Rcpp::export]]
Eigen::MatrixXd rInvWishRevCholesky_test(int v, Eigen::MatrixXd Psi){
  return rInvWishRevCholesky(v, Psi);
}

// [[Rcpp::export]]
Eigen::MatrixXd rMatUnitNormal_test1(int n, int m){
  MatrixXd X(n,m);
  fillUnitNormal(X);
  return X;
}

// [[Rcpp::export]]
Eigen::MatrixXd rMatUnitNormal_test2(int n){
  VectorXd X(n);
  fillUnitNormal(X);
  return X;
}
