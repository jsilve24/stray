context("test-pibble.R")
library(driver)
set.seed(4)


sim <- pibble_sim(true_priors=TRUE)


test_that("optim and uncollapse correctnesss", {
 
  init <- random_pibble_init(sim$Y)
  fit <- optimPibbleCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$KInv, 
                               sim$AInv, init,
                               n_samples=3000,
                               calcGradHess = FALSE, seed=40)
  
  # check closeness of MAP
  expect_true(abs(mean(fit$Pars - sim$Eta)) < .2)
  
  # Laplace approximation contains true value # given the true value
  # p0.25 <- apply(fit$Samples, c(1,2), function(x) quantile(x, probs=0.0025))
  # p99.75 <- apply(fit$Samples, c(1,2), function(x) quantile(x, probs=0.9975))
  #expect_true(sum(!((p0.25 <= Eta) & (p99.75 >= Eta))) < 0.2*N*(D-1))
  
  # Now check uncollapsing for Lambda
  fit2 <- uncollapsePibble(fit$Samples, sim$X, sim$Theta, sim$Gamma, 
                                     sim$Xi, sim$upsilon, seed=203)
  
  expect_true(mean(abs(apply(fit2$Lambda, c(1,2), mean) - sim$Phi)) < 0.7)
  p0.25 <- apply(fit2$Lambda, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit2$Lambda, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= sim$Phi) & (p99.75 >= sim$Phi))) < 0.01*sim$N*(sim$D-1))
  
  # check uncollapsing for Sigma
  # -- not implemented yet - correct results for Lambda imply correct results
  # -- for Sigma due to dependency structure. 
  
})


test_that("optim sylvester gets same result", {
  sim <- pibble_sim(D=30, N=10, true_priors = TRUE)
  init <- random_pibble_init(sim$Y)
  start <- Sys.time()
  fit <- optimPibbleCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$KInv, 
                               sim$AInv, init,
                               n_samples=2000,
                               calcGradHess = FALSE, 
                               useSylv=FALSE)
  end <- Sys.time()
  plain <- end-start
  start <- Sys.time()
  fitsylv <- optimPibbleCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$KInv, 
                               sim$AInv, init,
                               n_samples=2000,
                               calcGradHess = FALSE, 
                               useSylv=TRUE)
  end <- Sys.time()
  sylv <- end-start
  
  expect_equal(plain>sylv, TRUE)
  
  # check closeness of MAP
  expect_true(abs(mean(fit$Pars - fitsylv$Pars)) < .1)
})


test_that("pibble wrapper correctness", {
  fit <- pibble(sim$Y, sim$X, upsilon = sim$upsilon, Theta = sim$Theta, Xi=sim$Xi, 
                 Gamma=sim$Gamma, n_samples=3000)
  
  # Laplace approximation contains true value # given the true value
  p0.25 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= sim$Eta) & (p99.75 >= sim$Eta))) < 0.1*sim$N*(sim$D-1))
  
  # Check Lambda
  expect_true(mean(abs(apply(fit$Lambda, c(1,2), mean) - sim$Phi)) < 0.7)
  p0.25 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= sim$Phi) & (p99.75 >= sim$Phi))) < 0.05*sim$N*(sim$D-1))
  
  # check uncollapsing for Sigma
  # -- not implemented yet - correct results for Lambda imply correct results
  # -- for Sigma due to dependency structure. 
})


#' @param eta as an (D-1 x N x iter) array
uncollapse <- function(eta, X, upsilon, Theta, Xi, Gamma){
  d <- dim(eta)
  iter <- as.integer(d[3])
  N <- as.integer(d[2])
  D <- as.integer(d[1] + 1)
  Q <- as.integer(nrow(Gamma))
  Lambda <- array(0, c(D-1, Q, iter))
  Sigma <- array(0, c(D-1, D-1, iter))
  
  upsilonN <- upsilon+N
  GammaInv <- solve(Gamma)
  GammaN <- solve(tcrossprod(X)+ GammaInv)
  for (i in 1:iter){
    LambdaN <- (eta[,,i] %*% t(X) + Theta %*% GammaInv) %*% GammaN
    EN <- eta[,,i] - LambdaN %*% X
    Delta <- LambdaN - Theta
    XiN <- Xi + tcrossprod(EN) + Delta %*% solve(Gamma) %*% t(Delta)
    Sigma[,,i] <- MCMCpack::riwish(upsilonN, XiN)#solve(rWishart(1, upsilonN, XiN)[,,1])
    Z <- matrix(rnorm((D-1)*Q), D-1, Q)
    Lambda[,,i] <- LambdaN + t(chol(Sigma[,,i]))%*%Z%*%chol(GammaN)
  }
  return(list(Lambda=Lambda, Sigma=Sigma))
}

#' @param eta as an (D-1 x N x iter) array
uncollapse_mean_only <- function(eta, X, upsilon, Theta, Xi, Gamma){
  d <- dim(eta)
  iter <- as.integer(d[3])
  N <- as.integer(d[2])
  D <- as.integer(d[1] + 1)
  Q <- as.integer(nrow(Gamma))
  Lambda <- array(0, c(D-1, Q, iter))
  Sigma <- array(0, c(D-1, D-1, iter))
  
  upsilonN <- upsilon+N
  GammaInv <- solve(Gamma)
  GammaN <- solve(tcrossprod(X)+ GammaInv)
  for (i in 1:iter){
    LambdaN <- (eta[,,i] %*% t(X) + Theta %*% GammaInv) %*% GammaN
    Delta <- LambdaN - Theta
    EN <- eta[,,i] - LambdaN %*% X
    XiN <- Xi + tcrossprod(EN) + Delta %*% solve(Gamma) %*% t(Delta)
    Sigma[,,i] <- XiN*(upsilonN-D)
    Lambda[,,i] <- LambdaN 
  }
  return(list(Lambda=Lambda, Sigma=Sigma))
}


test_that("uncollapse correctnesss against double programming", {
  init <- random_pibble_init(sim$Y)
  fit <- optimPibbleCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$KInv, 
                               sim$AInv, init,
                               n_samples=2000,
                               calcGradHess = FALSE)
  
  # Now check uncollapsing for Lambda
  fit2 <- uncollapsePibble(fit$Samples, sim$X, sim$Theta, sim$Gamma, 
                                     sim$Xi, sim$upsilon, ret_mean = TRUE, 2234)
  
  dpres <- uncollapse_mean_only(fit$Samples, sim$X, sim$upsilon, sim$Theta, sim$Xi, 
                      sim$Gamma)
  
  expect_equal(fit2$Lambda, dpres$Lambda)
  expect_equal(fit2$Sigma, dpres$Sigma)
})


test_that("eigen and cholesky get same result", {
  sim <- pibble_sim(true_priors=TRUE, N=2, D=4)
  init <- random_pibble_init(sim$Y)
  fitc <- optimPibbleCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$KInv, 
                               sim$AInv, init,
                               n_samples=500000,
                               calcGradHess = FALSE, 
                               decomp="cholesky")
  fite <- optimPibbleCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$KInv, 
                                sim$AInv, init,
                                n_samples=500000,
                                calcGradHess = FALSE, 
                                decomp="eigen")
  
  expect_equal(apply(fitc$Samples, c(1,2), mean), 
               apply(fite$Samples, c(1,2), mean), 
               tolerance=0.01)
  
  expect_equal(apply(fitc$Samples, c(1,2), var), 
               apply(fite$Samples, c(1,2), var), 
               tolerance=0.01)
  
  expect_equal(fitc$logInvNegHessDet, fite$logInvNegHessDet, tolerance=1e-3)
})

test_that("logInvNegHessDet correct", {
  init <- random_pibble_init(sim$Y)
  fitc <- optimPibbleCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$KInv, 
                                sim$AInv, init,
                                n_samples=1,
                                calcGradHess = TRUE, 
                                decomp="cholesky")
  expect_equal(-log(det(-fitc$Hessian)), fitc$logInvNegHessDet)
})


test_that("max_iter leads to warning not error", {
  expect_warning(pibble(sim$Y, sim$X, max_iter=3))
})

