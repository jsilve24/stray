library(driver)
library(tidyverse)

D <- 10
Q <- 2
N <- 60

# Simulate Data
Sigma <- diag(sample(1:8, D-1, replace=TRUE))
Sigma[2, 3] <- Sigma[3,2] <- -1
Gamma <- matrix(c(1, .4, .4, 1.3), nrow=Q, byrow=TRUE)
Theta <- matrix(0, D-1, Q)
Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1)), nrow=D-1)%*%chol(Gamma)

X <- rbind(1, rnorm(N))
Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1)), nrow=D-1)
Pi <- t(alrInv(t(Eta)))
Y <- matrix(0, D, N)
for (i in 1:N) Y[,i] <- rmultinom(1, sample(500:10000), prob = Pi[,i])

# Priors
Xi <- diag(D-1)
upsilon <- D


K <- solve(Xi)
A <- solve(diag(N)+ t(X)%*%Gamma%*%X)

# fit <- optimMMTC(Y, upsilon, Theta%*%X, K, A, t(alr(t(Y)+0.65)), iter=2000)
# 



# foo <- chol(-fit$Hessian, pivot = TRUE)
# # 
# LH <- t(chol(solve(-fit$Hessian)))
# 
# z <- LH%*%matrix(rnorm(2000*30), 30, 2000)
# 
# 
# 




# double check hesssian
# 
# # Known working implementation of TVEC matrix
# TVEC <- function(m,n){
#   K <- matrix(0, m*n, m*n)
#   for (i in 1:(m*n)){
#     for(j in 1:(m*n)){
#       if (j == 1+m*(i-1)-(m*n-1)*floor((i-1)/n))  K[i,j] <- 1
#     }
#   }
#   return(K)
# }
# 
# 
# 
# E <- (fit$Pars - Theta%*%X)
# S <- diag(D-1) + K%*%E%*%A%*%t(E)
# O <- exp(fit$Pars)
# m <- 1+colSums(O)
# n <- colSums(Y)
# rho <- c(O)/ c(matrix(1, D-1, 1)%*%matrix(m, 1, N))
# R <- solve(S)%*%K
# C <- A%*%t(E)
# delta <- -0.5*(upsilon+N-D-2)
# grad <- t(c(Y[-D,])- c(matrix(1, D-1, 1)%*%n)*rho) # for multinomial
# grad <- grad + delta*t(c((R + t(R))%*%t(C)))
# grad[1:5]
# fit$Gradient[1:5]
# 
# L <- (C%*%R%*%t(C))%x%t(R)
# H <- matrix(0, N*(D-1), N*(D-1))
# for (j in 1:N){
#   it <- ((j-1)*(D-1)+1):(j*(D-1))
#   H[it, it] <- n[j]*(rho[it]%*%t(rho[it]) - diag(rho[it]))
# }
# tmp <- (A%x%(R+t(R))) - (L + t(L))
# # tmp <- tmp - TvecLMult(N, D-1,((R%*%t(C))%x%(C%*%t(R)) + (t(R)%*%t(C))%x%(C%*%R)))
# tmp <- tmp - TVEC(N, D-1)%*%((R%*%t(C))%x%(C%*%t(R)) + (t(R)%*%t(C))%x%(C%*%R))
# H <- H + delta*tmp
# fit$Hessian[1:5,1:5]
# H[1:5,1:5]
# 
# #eigen(solve(-H))
# LH <- t(chol(solve(-H)))
# z <- LH%*%matrix(rnorm(2000*30), 30, 2000)
