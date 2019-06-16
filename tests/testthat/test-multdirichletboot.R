context("test-multdirichletboot")
library(driver)

test_that("alr_default and alrInv_default is correct", {
  x <- miniclo_array(matrix(1:10, 5, 2), parts=1)
  x.alr <- alr_default_test(x)
  expect_equal(alr_array(x, parts=1), x.alr)
  expect_equal(x, alrInv_default_test(x.alr))
})

rDirichlet <- function(n_samples, alpha){
  p <- length(alpha)
  x <- matrix(0, p, n_samples)
  for (i in 1:p){
    x[i,] <- rgamma(n_samples, alpha[i], 1) 
  }
  return( miniclo_array(x, parts=1) )
}

test_that("MultDirichletBoot is correct", {
  n_samples <- 50000
  pi <- miniclo_array(matrix(1:5, 5, 1), parts=1)
  eta <- alr_array(pi, parts=1)
  depth <- 10
  Y <- matrix(rep(depth, 5), 5, 1)
  s <- MultDirichletBoot_test(n_samples, eta, Y, 0.05)
  
  x <- rDirichlet(n_samples, pi*depth*5)
  x <- alr_array(x, parts=1)
  
  expect_equal(rowMeans(x), rowMeans(s), tolerance=0.01)
  expect_equal(apply(x, 1, var), apply(s, 1, var), tolerance=0.05)
})

test_that("Timer does not have Error Johannes pointed out",{
  sim <- pibble_sim()
  fit <- pibble(sim$Y, sim$X, calcGradHess=FALSE, multDirichletBoot=0.65)
  expect_true(TRUE)
})
