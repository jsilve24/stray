context("test-matdist.R")

test_that("MVN handles non-square output", {
  Sigma <- matrix(c(1, -.5, -.5, 2),ncol=2, byrow = TRUE)
  E <- rMatNormalCholesky_test(matrix(0, 10, 2), diag(10), Sigma, discard=1)
  expect_equal(dim(E), c(10, 2))
})

test_that("MVN correctness of Mean", {
  M <- matrix(c(1,2,3,4), ncol=2)
  X2 <- array(0, dim=c(2, 2, 1000))
  for (i in 1:1000){
    X2[,,i] <- rMatNormalCholesky_test(M, diag(2), diag(2), discard=i)
  }

  # Tol based on Standard error of the mean
  expect_equal(apply(X2, c(1,2), mean), M, tolerance = 1/sqrt(1000)) 
})

test_that("MVN correctness of Covariances", {
  M <- matrix(0, ncol=2, nrow=2)
  U <- matrix(c(1,.5,.5, 2), ncol=2)
  V <- matrix(c(5,-.5,-.5, 1), ncol=2)
  LU <- t(chol(U))
  LV <- t(chol(V))
  t <- 10000
  X2 <- array(0, dim=c(2, 2, t))
  for (i in 1:t){
    X2[,,i] <- rMatNormalCholesky_test(M, LU,LV, discard=i)
  }
  
  # Test U
  tr <- function(x) sum(diag(x))
  
  # From Wikipedia
  # E[(X-M)(X-M)'] = Utr(V)
  X2_store <- array(0, dim=dim(X2))
  for (i in 1:t){
    X2_store[,,i] <- X2[,,i]%*%t(X2[,,i])
  }
  expect_equal(apply(X2_store, c(1,2), mean), U*tr(V), tolerance = 0.1)
  
  # E[(X-M)'(X-M)'] = Vtr(U)
  X2_store <- array(0, dim=dim(X2))
  for (i in 1:t){
    X2_store[,,i] <- t(X2[,,i])%*%(X2[,,i])
  }
  expect_equal(apply(X2_store, c(1,2), mean), V*tr(U), tolerance = 0.1)
})

test_that("InvWishart Correctness of Mean", {
  Psi <- matrix(c(1,.5,.5, 2), ncol=2)
  v <- 4
  t <- 100000
  Sigma <- array(0, dim=c(2,2,t))
  for (i in 1:t){
    Sigma[,,i] <- rInvWishRevCholesky_test(v, Psi)
    Sigma[,,i] <- tcrossprod(Sigma[,,i])
  }
  expect_equal(apply(Sigma, c(1,2), mean), Psi/(v-2-1), tolerance=0.1)
})


test_that("Unit normal filler is correct",{
  x <- rMatUnitNormal_test1(1000,1000)
  expect_equal(dim(x), c(1000, 1000))
  x <- c(x)
  expect_equal(mean(x), 0, tolerance=0.02)
  expect_equal(var(x), 1, tolerance=0.02)
  
  # Test that unit normal filler handles VectorXd objects as well. 
  x <- rMatUnitNormal_test2(100000)
  expect_equal(dim(x), c(100000, 1))
  x <- c(x)
  expect_equal(mean(x), 0, tolerance=0.02)
  expect_equal(var(x), 1, tolerance=0.02)
})

