context("test-iqlr")

test_that("iqlr gives correct format of output", {
  sim <- pibble_sim()
  fit <- pibble(sim$Y, sim$X)
  out <- lambda_to_iqlr(fit, 1:2)
  expect_equal(dim(out), c(ncategories(fit), ncovariates(fit), niter(fit)))
})

test_that("iqlr handles iter equal 1 case", {
  sim <- pibble_sim()
  fit <- pibble(sim$Y, sim$X,n_samples=1)
  fit <- pibblefit(as.integer(sim$D), as.integer(sim$N), as.integer(sim$Q), iter=1L,
                   coord_system="alr", alr_base=sim$D, 
                   Lambda = fit$Lambda)
  out <- lambda_to_iqlr(fit) # also testing handing of default focus.cov
  expect_equal(dim(out), c(ncategories(fit), ncovariates(fit), niter(fit)))
})