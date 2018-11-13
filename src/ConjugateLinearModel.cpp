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


//' Solve Bayesian Multivariate Conjugate Linear Model
//' 
//' See details for model.  Notation: \code{N} is number of samples,
//' \code{D} is the dimension of the response, \code{Q} is number
//' of covariates. 
//' 
//' @param Y matrix of dimension D x N
//' @param X matrix of covariates of dimension Q x N
//' @param Theta matrix of prior mean of dimension D x Q
//' @param Gamma covariance matrix of dimension Q x Q
//' @param Xi covariance matrix of dimension D x D
//' @param upsilon scalar (must be > D-1) degrees of freedom for InvWishart prior
//' @param n_samples number of samples to draw (default: 2000)
//' 
//' @details 
//'    \deqn{Y ~ MN_{D-1 x N}(Lambda*X, Sigma, I_N)}
//'    \deqn{Lambda ~ MN_{D-1 x Q}(Theta, Sigma, Gamma)}
//'    \deqn{Sigma ~ InvWish(upsilon, Xi)}
//' This function provides a means of sampling from the posterior distribution of 
//' \code{Lambda} and \code{Sigma}. 
//' @return List with components 
//' 1. Lambda Array of dimension (D-1) x Q x n_samples (posterior samples)
//' 2. Sigma Array of dimension (D-1) x (D-1) x n_samples (posterior samples)
//' @export
//' @md
//' @examples
//' sim <- mongrel_sim()
//' eta.hat <- t(driver::alr(t(sim$Y+0.65)))
//' fit <- conjugateLinearModel(eta.hat, sim$X, sim$Theta, sim$Gamma, 
//'                             sim$Xi, sim$upsilon, n_samples=2000)
// [[Rcpp::export]]
List conjugateLinearModel(const Eigen::Map<Eigen::MatrixXd> Y, 
                                const Eigen::Map<Eigen::MatrixXd> X, 
                                const Eigen::Map<Eigen::MatrixXd> Theta,
                                const Eigen::Map<Eigen::MatrixXd> Gamma, 
                                const Eigen::Map<Eigen::MatrixXd> Xi, 
                                const double upsilon, 
                                int n_samples = 2000){
  List out(2);
  out.names() = CharacterVector::create("Lambda", "Sigma");
  int Q = Gamma.rows();
  int D = Xi.rows();
  int N = X.cols();
  int iter = n_samples; // assumes result is an integer !!!
  double upsilonN = upsilon + N;
  MatrixXd GammaInv = Gamma.lu().inverse();
  MatrixXd GammaInvN = GammaInv + X*X.transpose(); 
  MatrixXd GammaN = GammaInvN.lu().inverse();
  MatrixXd LGammaN= GammaN.llt().matrixL();
  MatrixXd ThetaGammaInvGammaN = Theta*GammaInv*GammaN;
  MatrixXd XTGammaN = X.transpose()*GammaN;
  // // Storage for computation
  MatrixXd LambdaN(D, Q);
  MatrixXd XiN(D, D);
  MatrixXd LambdaDraw(D, Q);
  MatrixXd LSigmaDraw(D, D);
  MatrixXd SigmaDraw(D, D);
  MatrixXd ELambda(D, Q);
  MatrixXd EY(D, N);
  // Storage for output
  MatrixXd LambdaDrawO(D*Q, iter);
  MatrixXd SigmaDrawO(D*D, iter);

  // computation out of for-loop compared to mongrelcollapsed_uncollapse
  LambdaN = Y*XTGammaN+ThetaGammaInvGammaN;
  ELambda = LambdaN-Theta;
  EY = Y-LambdaN*X;
  XiN =  (EY*EY.transpose()).eval() + Xi + (ELambda*GammaInv*ELambda.transpose()).eval();
  
  // iterate over all draws of eta
  for (int i=0; i < iter; i++){
      R_CheckUserInterrupt();
      // Draw Random Component
      LSigmaDraw = rInvWishRevCholesky(upsilonN, XiN).matrix();
      // Note: correct even though LSigmaDraw is reverse cholesky factor
      LambdaDraw = rMatNormalCholesky(LambdaN, LSigmaDraw, LGammaN.matrix());

      // map output to vectors
      Map<VectorXd> LambdaDrawVec(LambdaDraw.data(), LambdaDraw.size());
      LambdaDrawO.col(i) = LambdaDrawVec;
      SigmaDraw = LSigmaDraw*LSigmaDraw.transpose();
      Map<VectorXd> SigmaDrawVec(SigmaDraw.data(), SigmaDraw.size());
      SigmaDrawO.col(i) = SigmaDrawVec;
  }

  IntegerVector dLambda = IntegerVector::create(D, Q, iter);
  IntegerVector dSigma = IntegerVector::create(D, D, iter);
  NumericVector nvLambda = wrap(LambdaDrawO);
  NumericVector nvSigma = wrap(SigmaDrawO);
  nvLambda.attr("dim") = dLambda;
  nvSigma.attr("dim") = dSigma;
  out[0] = nvLambda;
  out[1] = nvSigma;

  return out;
}
