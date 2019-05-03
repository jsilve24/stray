context("test-sylvester-speedups")

test_that("Sylvester Results Agree", {
  sim <- pibble_sim(D = 20, N=5)
  ThetaX <- sim$Theta %*% sim$X
  eta <- random_pibble_init(sim$Y)
  ll <- loglikPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
                               sylv=FALSE)
  llsylv <- loglikPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
                               sylv=TRUE)
  g <- gradPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
                            sylv=FALSE)
  gsylv <- gradPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
                                sylv=TRUE)
  hess <- hessPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
                            sylv=FALSE)
  hesssylv <- hessPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
                                sylv=TRUE)
  
  # microbenchmark::microbenchmark(gradPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
  #                                                     sylv=FALSE), 
  #                                gradPibbleCollapsed(sim$Y, sim$upsilon, ThetaX, sim$KInv, sim$AInv, eta, 
  #                                                     sylv=TRUE))
  expect_equal(ll, llsylv)
  expect_equal(g, gsylv)
  expect_equal(hess, hesssylv)
})
