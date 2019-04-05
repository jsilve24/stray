context("test-iqlr")

test_that("iqlr gives correct format of output", {
  sim <- pibble_sim()
  fit <- pibble(sim$Y, sim$X)
  out <- lambda_to_iqlr(fit, 1:2)
  expect_equal(dim(out), c(ncategories(fit), ncovariates(fit), niter(fit)))
})
