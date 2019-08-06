context("test-sylvester-speedups")

test_that("Pibble Sylvester Results Agree", {
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

test_that("Maltipoo Sylvester Results Agree", {
  # lazy, just use pibble sim data
  sim <- pibble_sim(D = 20, N=5)
  eta <- random_pibble_init(sim$Y)
  ell <- c(1)

  ll <- loglikMaltipooCollapsed(sim$Y, sim$upsilon, sim$Theta, sim$X, sim$KInv, sim$Gamma, eta, ell, 
                               sylv=FALSE)
  llsylv <- loglikMaltipooCollapsed(sim$Y, sim$upsilon, sim$Theta, sim$X, sim$KInv, sim$Gamma, eta, ell, 
                               sylv=FALSE)
  g <- gradMaltipooCollapsed(sim$Y, sim$upsilon, sim$Theta, sim$X, sim$KInv, sim$Gamma, eta, ell,
                            sylv=FALSE)
  gsylv <- gradMaltipooCollapsed(sim$Y, sim$upsilon, sim$Theta, sim$X, sim$KInv, sim$Gamma, eta, ell,
                            sylv=FALSE)
  hess <- hessMaltipooCollapsed(sim$Y, sim$upsilon, sim$Theta, sim$X, sim$KInv, sim$Gamma, eta, ell,
                            sylv=FALSE)
  hesssylv <- hessMaltipooCollapsed(sim$Y, sim$upsilon, sim$Theta, sim$X, sim$KInv, sim$Gamma, eta, ell,
                            sylv=FALSE)
  
  expect_equal(ll, llsylv)
  expect_equal(g, gsylv)
  expect_equal(hess, hesssylv)
})
