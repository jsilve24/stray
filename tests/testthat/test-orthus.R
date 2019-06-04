test_that("orthus sim and wrapper run without error", {
  sim <- orthus_sim()
  fit <- orthus(sim$Y, sim$Z, sim$X)
  expect_true(TRUE)
})

test_that("orthus wrapper correctness", {
  sim <- orthus_sim()
  fit <- orthus(sim$Y, sim$Z, sim$X, upsilon = sim$upsilon, Theta = sim$Theta, Xi=sim$Xi, 
                Gamma=sim$Gamma, n_samples=3000)
  
  # Laplace approximation contains true value # given the true value
  p0.25 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= sim$Eta) & (p99.75 >= sim$Eta))) < 0.1*sim$N*(sim$D-1))
  
  # Check Lambda
  expect_true(mean(abs(apply(fit$Lambda, c(1,2), mean) - sim$Phi)) < 0.5)
  p0.25 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= sim$Phi) & (p99.75 >= sim$Phi))) < 0.05*sim$N*(sim$D-1))
  
  # check uncollapsing for Sigma
  # -- not implemented yet - correct results for Lambda imply correct results
  # -- for Sigma due to dependency structure. 
})
