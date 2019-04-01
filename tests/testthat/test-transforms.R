context("test-transforms.R")


sim <- pibble_sim(D=4, Q=2, N=10, true_priors=TRUE)
fit <- pibble(sim$Y, sim$X)


test_that("pibble transform correctness", {
  ma <- to_alr(fit, 2)
  mc <- to_clr(fit)
  mi <- to_ilr(fit, driver::create_default_ilr_base(fit$D))
  expect_true(TRUE) # this is just here to get above to run
})
