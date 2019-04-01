context("test-laplaceapproximation")
library(driver)

test_that("eigen_lap gets correct answer", {
  n_samples <- 100000
  m <- 1:3
  S <- diag(4:6)
  S[1,2] <- S[2,1] <- -1
  #S <- -S
  z <- eigen_lap_test(n_samples, m, S, 0)
  
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
})


test_that("cholesky_lap gets correct answer", {
  n_samples <- 100000
  m <- 1:3
  S <- diag(4:6)
  S[1,2] <- S[2,1] <- -1
  #S <- -S
  z <- eigen_lap_test(n_samples, m, S, 0)
  
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
})

test_that("LaplaceApproximation gets correct result for full hessian", {
  n_samples <- 1000000
  m <- 1:3
  S <- diag(1:3)
  S[1,2] <- S[2,1] <- -1
  #S <- -S
  
  z <- LaplaceApproximation_test(n_samples, m, S, "eigen", 0)
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
  
  z <- LaplaceApproximation_test(n_samples, m, S, "cholesky", 0)
  expect_equal(var(t(z)), solve(S), tolerance=0.005)
  expect_equal(rowMeans(z), m, tolerance=.01)
})




