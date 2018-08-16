# library(MicrobeDS)
# library(phyloseq)
# #library(mongrel)
# library(dplyr)
# library(ape)
# 
# set.seed(91)
# 
# data("RISK_CCFA")
# 
# ## Huge
# dat <- RISK_CCFA %>%
#     subset_samples(disease_stat!="missing",
#                    immunosup!="missing") %>%
#     subset_samples(steroids=="false") %>%
#     subset_samples(antibiotics=="false") %>%
#     subset_samples(biologics=="false") %>%
#     prune_samples(sample_sums(.) >= 40000,.) %>%
#     filter_taxa(function(x) sum(x > 3) > 20, TRUE)
# 
# 
# sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>%
#   mutate(age = as.numeric(as.character(age)),
#          diseasesubtype=relevel(diseasesubtype, ref="no"),
#          disease_stat = relevel(disease_stat, ref="non-inflamed"))
# 
# X <- t(model.matrix(~disease_stat+age, data=sample_dat))
# Y <- otu_table(dat)
# 
# 
# upsilon <- ntaxa(dat)+3
# GG <- cbind(diag(ntaxa(dat)-1), -1)
# Xi <- (upsilon-ntaxa(dat)-2)*GG%*% diag(ntaxa(dat)) %*%t(GG) # update Xi prior
# Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
# Gamma <- diag(nrow(X))
# 
# 
# fit <- mongrel(Y, X, upsilon, Theta, Gamma, Xi, step_size=0.005, b1=0.99,
#                verbose=TRUE, verbose_rate=100, max_iter=200, jitter=1e-5,
#                decomp_method="eigen", n_samples=0)
