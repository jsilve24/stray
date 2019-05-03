# context("test-linesearch")
# 
# test_that("linesearch terminates with nonridiculous output", {
#   N <- 20
#   D <- 20
#   sim <- pibble_sim(D=D, N=N, true_priors=FALSE)
#   Z <- runif(N*(D-1))
#   ThetaX <- sim$Theta%*%sim$X
#   element <- 1
#   new_eta <- lineSearch(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, sim$Eta, element, 0.5, 0.0001)
#   vec_eta <- c(sim$Eta)
#   compare <- rep(TRUE, N*(D-1))
#   compare[element] <- FALSE
#   expect_equal(vec_eta[compare], new_eta[compare], tolerance=0.0001)
# })
# 
# 
# 
# 
