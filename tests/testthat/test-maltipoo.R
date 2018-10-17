context("test-maltipoo")

require(driver)
set.seed(4)

test_that("maltipoo wrapper correctness", {
  D <- 5; N <- 70; Q <- N
  
  X <- matrix(rnorm(N*Q), Q, N)
  delta_true <- .1
  U <- diag(Q)
  Gamma <- delta_true*U
  upsilon <- D+3000
  Xi <- diag(D-1)
  Sigma <- Xi/(upsilon-D)
  
  # Mean Zero
  Theta <- matrix(0, D-1, Q)
  Z <- matrix(rnorm(Q*(D-1)), D-1, Q)
  B <- Theta + t(chol(Sigma))%*%Z%*%chol(Gamma)
  Z <- matrix(rnorm(Q*(D-1)), D-1, N)
  Eta <- B%*%X + t(chol(Sigma))%*%Z
  Pi <- alrInv_array(Eta, coords=1)
  Y <- matrix(0, D, N)
  for (i in 1:N){
    Y[,i] <- rmultinom(1, 10000, prob=Pi[,i])
  }
  
  fit <- maltipoo(Y, X, upsilon, Theta, U, Xi, init=Eta, ellinit = log(delta_true))
  
  # Check that scale of VCs is correct
  expect_true(fit$VCScale-delta_true < 0.1)
  
  # Laplace approximation contains true value # given the true value
  p0.25 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= Eta) & (p99.75 >= Eta))) < 0.02*N*(D-1))
  
  # Check Lambda
  expect_true(mean(abs(apply(fit$Lambda, c(1,2), mean) - B)) < 0.01)
  p0.25 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= B) & (p99.75 >= B))) < 0.02*N*(D-1))
})
