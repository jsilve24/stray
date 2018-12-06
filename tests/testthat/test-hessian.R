context("test-hessian")

#' @param eta is  vec(eta) 
matt_nll <- function(eta, Y, X, upsilon, Theta, Xi, Gamma){
  D <- nrow(Y)
  N <- ncol(Y)
  eta <- matrix(eta, D-1, N)
  A <- solve(diag(N) + t(X)%*%Gamma %*% X)
  K <- solve(Xi)
  delta <- (upsilon+N+D-2)/2
  E <- eta - Theta%*%X
  S <- diag(D-1) + K%*%E%*%A%*%t(E)
  nll <- delta*log(det(S))
  return(nll)
}

#' @param eta is  vec(eta) 
mult_nll <- function(eta, Y, X, upsilon, Theta, Xi, Gamma){
  D <- nrow(Y)
  N <- ncol(Y)
  nll <- 0
  eta <- matrix(eta, D-1, N)
  pi <- t(driver::alrInv(t(eta)))
  for (i in 1:N){
    nll <- nll - dmultinom(Y[,i], prob=pi[,i], log=TRUE) # NEGATIVE!
  }
  return(nll)
}

#' @param eta is  vec(eta) 
nll <- function(eta, Y, X, upsilon, Theta, Xi, Gamma){
  nll <- 0
  nll <- nll+ matt_nll(eta, Y, X, upsilon, Theta, Xi, Gamma)
  nll <- nll+mult_nll(eta, Y, X, upsilon, Theta, Xi, Gamma)
  return(nll)
}

#' function to calculate hessian for model at a given eta 
#' @param eta (D-1) x N matrix 
#' @details uses hessMongrelCollapsed function of mongrel 
hessMC <- function(mdataset, eta){
  X <- mdataset$X
  A <- solve(diag(mdataset$N)+ t(X)%*%mdataset$Gamma%*%X)
  hessMongrelCollapsed(mdataset$Y, mdataset$upsilon, 
                       mdataset$Theta%*%X, solve(mdataset$Xi), 
                       A, eta)
}

sim <- mongrel_sim(D=5, N=10, true_priors=FALSE)
nll_partial <- function(x) nll(x, sim$Y, sim$X, sim$upsilon,
                               sim$Theta, sim$Xi, sim$Gamma)

test_that("hessian agrees with finite differences", {
  hess.nd <- numDeriv::hessian(nll_partial, c(sim$Eta))
  A <- solve(diag(sim$N) + t(sim$X) %*% sim$Gamma %*% sim$X)
  hess <- hessMongrelCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, 
                               solve(sim$Xi), A, sim$Eta)
  expect_equal(hess.nd, -hess, tolerance=1e-3)
})
