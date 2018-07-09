context("test-main.R")
library(driver)
library(tidyverse)


sim <- mongrel_sim(true_priors=TRUE)
attach(sim)


test_that("optim and uncollapse correctnesss", {
 
  init <- random_mongrel_init(Y)
  fit <- optimMongrelCollapsed(Y, upsilon, Theta%*%X, K, A, init,
                               n_samples=2000,
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
  fit <- fit_mongrel(Y, X)
  
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

