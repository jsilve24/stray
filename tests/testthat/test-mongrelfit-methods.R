context("test-pibblefit-methods.R")

sim <-   pibble_sim()

test_that("sample_prior correct", {
  fit <- pibble(sim$Y, sim$X)
  priors <- sample_prior(fit)
  
  expect_equal(apply(priors$Eta, c(1, 2), mean), priors$Theta%*%priors$X, 
               tolerance=0.1)
  
  expect_equal(apply(priors$Lambda, c(1, 2), mean), priors$Theta, 
               tolerance=0.1)
})


test_that("Predict works with priors only",{
  fit <- pibble(Y=NULL, sim$X,  D=sim$D)
  foo <- predict(fit, response="Y", size = 5000)
  expect_equal(dim(foo), c(sim$D, sim$N, fit$iter))
})

# Fixing github issue #4
test_that("Plot works with focus.cooord and coord system change",{
  fit <- pibble(sim$Y, sim$X)
  fit <- to_clr(fit)
  p <- plot(fit, par="Lambda", focus.coord=c("clr_c1",  "clr_c2"))
  #p <- plot(fit, par="Lambda", focus.coord=1:4) # Not yet implemented
  expect_error(print(p), NA)
})


test_that("Predict works in CLR", {
  fit <- pibble(sim$Y, sim$X)
  fit <- to_clr(fit)
  expect_error(predict(fit), NA) 
})

test_that("summary works with pars=Lambda", {
  fit <- pibble(sim$Y, sim$X, pars="Lambda")
  summary(fit)
  expect_true(TRUE)
})


# test_that("Plots work with iter=1", {
#   fit <- pibble(sim$Y, sim$X, n_samples=1)
#   plot(fit, par="Sigma")
#   plot(fit, par="Eta")
# })