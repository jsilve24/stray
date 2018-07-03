context("test-transforms.R")

D <- 4
Q <- 2
N <- 10

# Simulate Data
Sigma <- diag(sample(1:8, D-1, replace=TRUE))
Sigma[2, 3] <- Sigma[3,2] <- -1
Gamma <- diag(sqrt(rnorm(Q)^2))
Theta <- matrix(0, D-1, Q)
Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1)), nrow=D-1)%*%chol(Gamma)
X <- matrix(rnorm(N*(Q-1)), Q-1, N)
X <- rbind(1, X)
Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1)), nrow=D-1)
Pi <- t(alrInv(t(Eta)))
Y <- matrix(0, D, N)
for (i in 1:N) Y[,i] <- rmultinom(1, sample(5000:10000), prob = Pi[,i])

# Priors
#upsilon <- D
#Xi <- diag(D-1)
upsilon <- D+10
Xi <- Sigma*(upsilon-D-2)

# Precompute
K <- solve(Xi)
A <- solve(diag(N)+ t(X)%*%Gamma%*%X)
fit <- mongrel(Y, X)

test_that("mongrel transform correctness", {
  ma <- mongrel_to_alr(fit, 2)
  mc <- mongrel_to_clr(fit)
  mi <- mongrel_to_ilr(fit, driver::create_default_ilr_base(fit$D))
  
  # more to come... 
  
  expect_error(mongrel_to_alr(ma, 2))
})
