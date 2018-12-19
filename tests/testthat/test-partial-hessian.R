context("test-partial_hessian.R")

library(driver)

D <- 50
N <- 10
Q <- 2
sim <- mongrel_sim(N=N, D=D, Q=Q, true_priors=TRUE)

test_that("mongrel partial hessian optim runs", {
  init <- random_mongrel_init(sim$Y)
  
  # I want to remember that this is broken (gives shitty result and crashes 
  # with MKL Settings)
  expect_true(FALSE)
  # fit <- optimMongrelCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$K, 
  #                              sim$A, init,
  #                              n_samples=1000,
  #                              calcGradHess = TRUE,
  #                              decomp_method="eigen",
  #                              calcPartialHess = TRUE)

  # mult_t_hessian <- function(x, nu, m, K, omegainv){
  #   p <- nrow(x)
  #   e <- x-m
  #   delta <- (nu+p)/2
  #   s <- 1+omegainv*t(e)%*%K%*%e
  #   C <- K%*%e
  #   R <- omegainv/s
  #   L <- c(R)^2*tcrossprod(C)
  #   return(delta*(2*c(R)*K - 4*L))
  # }
  # x <- fit$Pars[,1,drop=F]
  # K <- sim$K
  # m <- (sim$Theta%*%sim$X)[,1,drop=F]
  # nu <- sim$upsilon
  # omegainv <- sim$A[1,1]
  # 
  # foo <- mult_t_hessian(x,nu, m, K,omegainv)
  # -foo[1:5,1:5] - fit$Hessian[1:49,1:49][1:5,1:5]
  # 
  # eigen(fit$Hessian[1:49,1:49])
  
  expect_true(TRUE)
})

test_that("mongrel multDirichletBoot optim runs", {
  init <- random_mongrel_init(sim$Y)
  fit <- optimMongrelCollapsed(sim$Y, sim$upsilon, (sim$Theta%*%sim$X), sim$K, 
                               sim$A, init,
                               n_samples=2000,
                               multDirichletBoot = TRUE)

  sim$Eta[1:5,1:5]
  apply(fit$Samples, c(1,2), mean)[1:5,1:5]
  expect_true(TRUE)
})