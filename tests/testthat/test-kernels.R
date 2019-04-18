context("test-kernels")

test_that("RBF works", {
  X <- matrix(rnorm(10), 2, 5)
  G <- RBF(X, 2, .2)
  expect_true(all(eigen(G)$values>0))
  expect(all(diag(G)-4 < .000001))
})

test_that("LINEAR works", {
  X <- matrix(rnorm(15), 5, 3)
  G <- LINEAR(X, 1, rep(0, nrow(X)))
  expect_true(all(eigen(G)$values>0))
})