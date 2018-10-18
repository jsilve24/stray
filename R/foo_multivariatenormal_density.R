# N <- 100
# y <- rnorm(N, 0, 5)
# 
# ll <- function(sigma){
#   #-0.5*sum((1/sigma)*y^2) - 0.5*log(2*pi*N*sigma)
#   sum(dnorm(y, 0, sd = sqrt(sigma), log = TRUE))
# }
# plot(1:100, sapply(1:100,ll))
# 
# 
# 
# ll <- function(sigma){
#   Sigma <- sigma*diag(N)
#   -0.5*log(det(2*pi*Sigma)) - 0.5*t(y)%*%solve(Sigma)%*%y
# }
# plot(1:100, sapply(1:100,ll))
