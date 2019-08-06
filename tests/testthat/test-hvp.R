# context("test-hessianvectorproduct")
# 
# 
# test_that("hessVectorProd output is reasonable", {
#   set.seed(88)
#   N <- 20
#   D <- 20
#   sim <- pibble_sim(D=D, N=N, true_priors=FALSE)
#   Z <- runif(N*(D-1))
#   ThetaX <- sim$Theta%*%sim$X
#   prod1 <- hessPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, sim$Eta)
#   prod1 <- prod1%*%Z;
#   prod2 <- hessVectorProd(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, sim$Eta, Z, 0.001)
#   expect_equal(prod1[,1], prod2, tolerance=1e-5)
# })
# 
# 


