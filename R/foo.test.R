# f <- function(eta, Y){
#   pi <- alrInv(c(eta))
#   dmultinom(Y, prob=pi, log = TRUE)
# }
# 
# 
# library(mongrel)
# library(tidyverse)
# attach(sim_data)
# fit <- list()
# for (i in 1:100){
#   Yboot <- Y; Yboot[] <- 0
#   for (j in 1:sim_data$N){
#     Yboot[,j] <- rmultinom(1, size=sum(Y[,j]), c(miniclo(c(Y[,j]+.5))))
#   }
#   fit[[i]] <- mongrel(Y=Yboot, X=X, upsilon, Theta, Gamma, Xi, init=random_mongrel_init(Y), n_samples=0)  
# }
# 
# fit %>% 
#   map("Lambda") %>% 
#   map(gather_array, val, coordinate, covariate, dim_3) %>% 
#   bind_rows(.id="iter") %>% 
#   select(-dim_3) %>% 
#   group_by(coordinate, covariate) %>% 
#   summarise_posterior(val)
#   
# 
# 
# 
# 
# 
# 
# 
# f.partial <- function(x) f(x, Y[,3])
# foo <- numDeriv::hessian(f.partial,  fit$Eta[,3,1])
# foo <- solve(foo)
# 
# 
# # dirichlet symmetry? -----------------------------------------------------
# 
# rDirichlet <- function(n, alpha){
#   unclass(compositions::rDirichlet.acomp(n, alpha))
# }
# 
# foo <- rDirichlet(10000, Y[,3]+.)
# foo <- alr(foo)
# foo[, 21] %>% density() %>% plot()
