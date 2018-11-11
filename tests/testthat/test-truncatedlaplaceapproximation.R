context("test-truncatedlaplaceapproximation")

test_that("Truncated Eigen Decomposition Correct", {
  N <- 30
  D <- 30
  sim <- mongrel_sim(N=N, D=D, true_priors = TRUE)
  etavec <- rnorm(N*(D-1))
  eigs <- MongrelTruncatedEigen_mongrel_test(sim$Y, sim$upsilon, 
                                             sim$Theta%*%sim$X, 
                                             sim$K, sim$A, etavec, 
                                             0.001, 20, 100)
  eta <- matrix(etavec, D-1, N)
  hess <- hessMongrelCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$K, sim$A, eta)
  
  # microbenchmark::microbenchmark( eigs <- MongrelTruncatedEigen_mongrel_test(sim$Y, sim$upsilon, 
  #                                                                             sim$Theta%*%sim$X, 
  #                                                                             sim$K, sim$A, etavec, 
  #                                                                             0.001, 2, 30), 
  #                                 eigs <- MongrelTruncatedEigen_mongrel_test(sim$Y, sim$upsilon, 
  #                                                                            sim$Theta%*%sim$X, 
  #                                                                            sim$K, sim$A, etavec, 
  #                                                                            0.001, 1, 30),
  #                                 eigs <- MongrelTruncatedEigen_mongrel_test(sim$Y, sim$upsilon, 
  #                                                                            sim$Theta%*%sim$X, 
  #                                                                            sim$K, sim$A, etavec, 
  #                                                                            0.001, 10, 600),
  #                                 #{  hess <- hessMongrelCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$K, sim$A, eta);
  #                                 #foo <- eigen(-hess)}, 
  #                                 chol(-hess, pivot=TRUE), 
  #                                 times = 1)
  foo <- eigen(-hess)
  expect_equal(rev(foo$vals)[1:20], eigs$vals)

  
})
