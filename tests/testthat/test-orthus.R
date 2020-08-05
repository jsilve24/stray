set.seed(859)

test_that("orthus sim and wrapper run without error", {
  sim <- orthus_sim()
  fit <- orthus(sim$Y, sim$Z, sim$X)
  expect_true(TRUE)
})

test_that("orthus wrapper correctness", {
  sim <- orthus_sim()
  fit <- orthus(sim$Y, sim$Z, sim$X, upsilon = sim$upsilon, Theta = sim$Theta, Xi=sim$Xi, 
                Gamma=sim$Gamma, n_samples=3000)
  
  # Laplace approximation contains true value # given the true value
  p0.25 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Eta, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= sim$Eta) & (p99.75 >= sim$Eta))) < 0.2*sim$N*(sim$D-1))
  
  # Check Lambda
  expect_true(mean(abs(apply(fit$Lambda, c(1,2), mean) - sim$Phi)) < 0.5)
  p0.25 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.0025))
  p99.75 <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=0.9975))
  expect_true(sum(!((p0.25 <= sim$Phi) & (p99.75 >= sim$Phi))) < 0.05*sim$N*(sim$D-1))

})


test_that("Orthus works with multDirichletBoot", {
  sim <- orthus_sim()
  fit <- orthus(sim$Y, sim$Z, sim$X, upsilon = sim$upsilon, Theta = sim$Theta, Xi=sim$Xi, 
                Gamma=sim$Gamma, multDirichletBoot=.5)
  expect(TRUE, "cannot fail")
})

# test_that("orthus identical results with fixed seed", {
#   set.seed(3)
#   sim <- orthus_sim()
#   fit <- orthus(sim$Y, sim$Z, sim$X, seed=5)
#   Lambda.test <- fit$Lambda[1:5,1:2,1:5]
#   Lambda.test[1:5,,1]
#   #save(Lambda, file="tests/Lambda_seed3-5.RData")
#   load("tests/Lambda_seed3-5.RData")
#   expect_equal(Lambda.test, Lambda)
# })



