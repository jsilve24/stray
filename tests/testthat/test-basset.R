context("test-basset")

set.seed(569)

test_that("basset and predict.bassetfit run", {
  sim <- pibble_sim()
  Gamma <- function(X) SE(X)
  Theta <- function(X) matrix(0, nrow(sim$Y)-1, ncol(X))
  fit <- basset(sim$Y, sim$X, Gamma = Gamma, Theta = Theta)
  foo <- predict(fit, matrix(c(1,2)))
  expect_true(TRUE)
})

