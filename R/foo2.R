mult_t_density <- function(x, nu, m, K, omegainv){
  p <- nrow(x)
  e <- x-m
  delta <- (nu + p)/2
  delta*log(1+omegainv*t(e)%*%K%*%e)
}

mat_t_density <- function(X, nu, M, K, A){
  p <- nrow(X)
  N <- ncol(X)
  E <- X-M
  delta <- (nu+N + p -1)/2
  delta*log(det(diag(p) + K%*%E%*%A%*%t(E)))
}

nu <-10
p <- 4
N <- 1
M <- matrix(0, p, N)
K <- rWishart(1, 10, diag(p))[,,1]
#A <- diag(diag(rWishart(1, 10, diag(N))[,,1]))
A <- 4
X <- matrix(rnorm(p*N), p, N)

# Test 1 
foo <- rep(0, N)
for (i in 1:N){
  foo[i] <- mult_t_density(X[,i,drop=F], nu, M[,i,drop=F], K, diag(A)[i])
}
sum(foo)
mat_t_density(X, nu, M, K, A)


# Test 2 - check hessian
mult_t_hessian <- function(x, nu, m, K, omegainv){
  p <- nrow(x)
  e <- x-m
  delta <- (nu+p)/2
  s <- 1+omegainv*t(e)%*%K%*%e
  C <- K%*%e
  R <- omegainv/s
  L <- c(R)^2*tcrossprod(C)
  return(delta*(2*c(R)*K - 4*L))
}

numDeriv::hessian(function(x) mult_t_density(x, nu, M, K, A), X)
mult_t_hessian(X, nu, M, K, A)
