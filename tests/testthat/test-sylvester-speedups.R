context("test-sylvester-speedups")

test_that("Sylvester Results Agree", {
  sim <- mongrel_sim(D = 20, N=5)
  ThetaX <- sim$Theta %*% sim$X
  eta <- random_mongrel_init(sim$Y)
  ll <- loglikMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
                               sylv=FALSE)
  llsylv <- loglikMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
                               sylv=TRUE)
  g <- gradMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
                            sylv=FALSE)
  gsylv <- gradMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
                                sylv=TRUE)
  hess <- hessMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
                            sylv=FALSE)
  hesssylv <- hessMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
                                sylv=TRUE)
  
  # microbenchmark::microbenchmark(gradMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
  #                                                     sylv=FALSE), 
  #                                gradMongrelCollapsed(sim$Y, sim$upsilon, ThetaX, sim$K, sim$A, eta, 
  #                                                     sylv=TRUE))
  expect_equal(ll, llsylv)
  expect_equal(g, gsylv)
  expect_equal(hess, hesssylv)
})
