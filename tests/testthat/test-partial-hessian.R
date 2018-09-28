context("test-partial_hessian.R")

library(driver)
set.seed(4)

D <- 6
N <- 10
Q <- 2
sim <- mongrel_sim(N=N, D=D, Q=Q, true_priors=TRUE)

test_that("stub", {
  init <- random_mongrel_init(sim$Y)
  fit <- optimMongrelCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$K, 
                               sim$A, init,
                               n_samples=1000,
                               calcGradHess = TRUE,
                               decomp_method="eigen",
                               calcPartialHess = TRUE)

  # expect_equal(fit1[3], fit2[3], tolerance = 0.1)

  expect_true(TRUE)
})