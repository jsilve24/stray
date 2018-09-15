# #' Simulate Data and Priors for Maltipoo
# #' 
# #' @inheritParams mongrel_sim
# #' @param P number of variance components to simulate
# #' 
# #' @return list
# #' @export
# #' 
# #' @examples
# #' maltipoo_sim(10, 30, 30, 3, TRUE, FALSE)
# maltipoo_sim <- function(D=10, N=30, P=3,
#                          use_names=TRUE, true_priors=FALSE){
# 
#   Q <- N
# 
#   # Simulate Data
#   Sigma <- diag(sample(1:8, D-1, replace=TRUE))
#   Sigma[2, 3] <- Sigma[3,2] <- -1
# 
#   rU <- rWishart(P, Q+3, diag(Q))
#   U <- matrix(0, P*Q, Q)
#   Gamma_true <- matrix(0, Q, Q)
#   VCScale_true <- rgamma(P, 1, 1)
#   for (i in 1:P){
#     U[((i-1)*Q+1):(i*Q), ] <- solve(rU[,,i])
#     Gamma_true <- Gamma_true + VCScale_true[i]^2*U[((i-1)*Q+1):(i*Q), ]
#   }
#   rm(rU)
# 
#   Gamma_true <- diag(sqrt(rnorm(Q)^2))
#   Theta <- matrix(0, D-1, Q)
#   Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1)), nrow=D-1)%*%chol(Gamma_true)
#   X <- diag(Q)
#   #X <- rbind(1, X)
#   Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1)), nrow=D-1)
#   Pi <- t(driver::alrInv(t(Eta)))
#   Y <- matrix(0, D, N)
#   for (i in 1:N) Y[,i] <- rmultinom(1, sample(5000:10000), prob = Pi[,i])
#   if (use_names){
#     colnames(X) <- colnames(Y) <- paste0("s", 1:N)
#     rownames(Y) <- paste0("c", 1:D)
#     rownames(X) <- paste0("x", 1:Q)
#   }
# 
#   # Priors
#   if (true_priors){
#     upsilon <- D+50
#     Xi <- Sigma*(upsilon-D-2)
#   } else {
#     upsilon <- D
#     Xi <- diag(D-1)
#   }
# 
#   # Precompute
#   K <- solve(Xi)
# 
#   return(list(Sigma=Sigma, Gamma_true=Gamma_true, D=D, N=N, Q=Q, Theta=Theta, Phi=Phi,
#               X=X, Y=Y, Eta=Eta, upsilon=upsilon, Xi=Xi, K=K, U=U, VCScale_true=VCScale_true))
# }
