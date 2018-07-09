context("test-transforms.R")


sim <- mongrel_sim(D=4, Q=2, N=10, true_priors=TRUE)
fit <- mongrel(sim$Y, sim$X)


test_that("mongrel transform correctness", {
  ma <- mongrel_to_alr(fit, 2)
  mc <- mongrel_to_clr(fit)
  mi <- mongrel_to_ilr(fit, driver::create_default_ilr_base(fit$D))
  expect_true(TRUE) # this is just here to get above to run
})
