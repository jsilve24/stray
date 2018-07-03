context("test-main.R")
library(driver)
library(tidyverse)

D <- 10
Q <- 2
N <- 30

# Simulate Data
Sigma <- diag(sample(1:8, D-1, replace=TRUE))
Sigma[2, 3] <- Sigma[3,2] <- -1
Gamma <- diag(sqrt(rnorm(Q)^2))
Theta <- matrix(0, D-1, Q)
Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1)), nrow=D-1)%*%chol(Gamma)
X <- matrix(rnorm(N*(Q-1)), Q-1, N)
X <- rbind(1, X)
Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1)), nrow=D-1)
Pi <- t(alrInv(t(Eta)))
Y <- matrix(0, D, N)
for (i in 1:N) Y[,i] <- rmultinom(1, sample(5000:10000), prob = Pi[,i])
colnames(X) <- colnames(Y) <- paste0("s", 1:N)
rownames(Y) <- paste0("c", 1:D)
rownames(X) <- paste0("x", 1:Q)


# Priors
#upsilon <- D
#Xi <- diag(D-1)
upsilon <- D+10
Xi <- Sigma*(upsilon-D-2)

# Precompute
K <- solve(Xi)
A <- solve(diag(N)+ t(X)%*%Gamma%*%X)


test_that("optim and uncollapse correctnesss", {
 
  init <- random_mongrel_init(Y)
  fit <- optimMongrelCollapsed(Y, upsilon, Theta%*%X, K, A, init,
                               iter=2000, numexcessthresh=0,
                               calcGradHess = FALSE)
  
  # check closeness of MAP
  expect_true(abs(mean(fit$Pars - Eta)) < .5)
  
  # Laplace approximation contains true value # given the true value
  p0.25 <- apply(fit$Samples, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Samples, c(1,2), function(x) quantile(x, probs=0.9975))
  #expect_true(sum(!((p0.25 <= Eta) & (p99.75 >= Eta))) < 0.2*N*(D-1))
  
  # Now check uncollapsing for Lambda
  fit2 <- uncollapseMongrelCollapsed(fit$Samples, X, Theta, Gamma, Xi, upsilon)
  
  expect_true(mean(abs(apply(fit2$Lambda, c(1,2), mean) - Phi)) < 0.5)
  p0.25 <- apply(fit2$Lambda, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit2$Lambda, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= Phi) & (p99.75 >= Phi))) < 0.2*N*(D-1))
  
  # check uncollapsing for Sigma
  # -- not implemented yet - correct results for Lambda imply correct results
  # -- for Sigma due to dependency structure. 
  
})


test_that("mongrel wrapper correctness", {
  fit <- mongrel(Y, X)
  
  # Laplace approximation contains true value # given the true value
  p0.25 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= Eta) & (p99.75 >= Eta))) < 0.2*N*(D-1))
  
  # Check Lambda
  expect_true(mean(abs(apply(fit$Lambda, c(1,2), mean) - Phi)) < 0.5)
  p0.25 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= Phi) & (p99.75 >= Phi))) < 0.2*N*(D-1))
  
  # check uncollapsing for Sigma
  # -- not implemented yet - correct results for Lambda imply correct results
  # -- for Sigma due to dependency structure. 
})
