context("test-sample_prior-correct.R")

test_that("sample_prior correct", {
  sim <- mongrel_sim()
  fit <- mongrel(sim$Y, sim$X)
  priors <- sample_prior(fit)
  
  expect_equal(apply(priors$Eta, c(1, 2), mean), priors$Theta%*%priors$X, 
               tolerance=0.1)
  
  expect_equal(apply(priors$Lambda, c(1, 2), mean), priors$Theta, 
               tolerance=0.1)
})
