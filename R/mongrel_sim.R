#' Simulate simple mongrel dataset and priors (for testing)
#' 
#' @param D number of multinomial categories
#' @param N number of samples
#' @param Q number of covariates (first one is an intercept, must be > 1)
#' @param use_names should samples, covariates, and categories be named
#' @param true_priors should Xi and upsilon be chosen to have mean at true 
#'   simulated value
#' @return list
#' @export
#' @importFrom driver alrInv
#' @importFrom stats rnorm rmultinom
#' @examples 
#' sim <- mongrel_sim()
mongrel_sim <- function(D=10, N=30, Q=2, use_names=TRUE, true_priors=FALSE){
  
  # Simulate Data
  Sigma <- diag(sample(1:8, D-1, replace=TRUE))
  Sigma[2, 3] <- Sigma[3,2] <- -1
  Gamma <- diag(sqrt(rnorm(Q)^2))
  Theta <- matrix(0, D-1, Q)
  Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1)), nrow=D-1)%*%chol(Gamma)
  X <- matrix(rnorm(N*(Q-1)), Q-1, N)
  X <- rbind(1, X)
  Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1)), nrow=D-1)
  Pi <- t(driver::alrInv(t(Eta)))
  Y <- matrix(0, D, N)
  for (i in 1:N) Y[,i] <- rmultinom(1, sample(5000:10000), prob = Pi[,i])
  if (use_names){
    colnames(X) <- colnames(Y) <- paste0("s", 1:N)
    rownames(Y) <- paste0("c", 1:D)
    rownames(X) <- paste0("x", 1:Q)  
  }
  
  # Priors
  if (true_priors){
    upsilon <- D+10
    Xi <- Sigma*(upsilon-D)
  } else {
    upsilon <- D
    Xi <- diag(D-1)
  }
  
  # Precompute
  K <- solve(Xi)
  A <- solve(diag(N)+ t(X)%*%Gamma%*%X)
  
  return(list(Sigma=Sigma, Gamma=Gamma, D=D, N=N, Q=Q, Theta=Theta, Phi=Phi, 
              X=X, Y=Y, Eta=Eta, upsilon=upsilon, Xi=Xi, K=K, A=A))
}