# N <- 20; D <- 40
# upsilon <- D-1+10
# delta <- 3
# Sigma <- solve(rWishart(1, upsilon, diag(D-1))[,,1])
# Eta <- matrix(0, D-1, N)
# for (i in 1:N){  
#   Z <- matrix(rnorm(D-1), D-1, 1)
#   Eta[,i] <- delta*t(chol(Sigma))%*%Z
# }
# 
# 
# ll <- function(delta){
#   Omega <- delta^2*diag(N)
#   S <- diag(D-1) + solve(diag(D-1))%*%(Eta)%*%solve(Omega)%*%t(Eta)
#   -(D-1)/2*log(det(Omega)) - 0.5*(upsilon+N+D-2)*log(det(S))
# }
# 
# plot(1:20, sapply(1:20, ll))
# 
# 
# # with regression ---------------------------------------------------------
# 
# D <- 5; N <- 70; Q <- N
# 
# X <- matrix(rnorm(N*Q), Q, N)
# delta_true <- 5
# U <- diag(Q)
# Gamma <- delta_true*U
# upsilon <- D+30
# Xi <- diag(D-1)
# Sigma <- Xi/(upsilon-D-2)
# 
# # Mean Zero
# Theta <- matrix(0, D-1, Q)
# Z <- matrix(rnorm(Q*(D-1)), D-1, Q)
# B <- Theta + t(chol(Sigma))%*%Z%*%chol(Gamma)
# Z <- matrix(rnorm(Q*(D-1)), D-1, N)
# Eta <- B%*%X + t(chol(Sigma))%*%Z
# 
# 
# ll <- function(delta){
#   E <- Eta-Theta%*%X
#   K <- solve(Xi)
#   XTUX <- t(X)%*%U%*%X
#   Ainv <- diag(N) + delta*XTUX
#   A <- solve(Ainv)
#   S <- diag(D-1)+K%*%E%*%A%*%t(E)
#   const <- 0.5*(upsilon+N+D-2) 
#   ll <- -const*log(det(S)) - (0.5*(D-1))*log(det(Ainv))
#   return(ll)
# }
# 
# plot(1:20, sapply(1:20, ll))
# 
# which.max(sapply(1:20, ll))
# 
# 
# # with exponential parameterization regression -------------------------------
# 
# D <- 5; N <- 70; Q <- N
# 
# X <- matrix(rnorm(N*Q), Q, N)
# delta_true <- 5
# U <- diag(Q)
# Gamma <- delta_true*U
# upsilon <- D+30
# Xi <- diag(D-1)
# Sigma <- Xi/(upsilon-D-2)
# 
# # Mean Zero
# Theta <- matrix(0, D-1, Q)
# Z <- matrix(rnorm(Q*(D-1)), D-1, Q)
# B <- Theta + t(chol(Sigma))%*%Z%*%chol(Gamma)
# Z <- matrix(rnorm(Q*(D-1)), D-1, N)
# Eta <- B%*%X + t(chol(Sigma))%*%Z
# 
# 
# ll <- function(ell){
#   E <- Eta-Theta%*%X
#   K <- solve(Xi)
#   XTUX <- t(X)%*%U%*%X
#   Ainv <- diag(N) + exp(ell)*XTUX
#   A <- solve(Ainv)
#   S <- diag(D-1)+K%*%E%*%A%*%t(E)
#   const <- 0.5*(upsilon+N+D-2) 
#   ll <- -const*log(det(S)) - (0.5*(D-1))*log(det(Ainv))
#   return(ll)
# }
# 
# plot(log(1:20), sapply(log(1:20), ll))
# 
# which.max(sapply(log(1:20), ll))
