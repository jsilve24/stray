#' Transform Lambda into IQLR (Inter-Quantile Log-Ratio)
#' 
#' Primarily intended for doing differential expression analysis under 
#' assumption that only small group of categories (e.g., taxa / genes) are changing 
#' 
#' @param m object of class pibblefit (e.g., output of \code{\link{pibble}})
#' @param focus.cov vector of integers or characters specifying columns (covariates)
#'   of Lambda to include in calculating IQLR (if NULL, default, then uses all covariates)
#' @param probs bounds for categories (i.e., features / genes / taxa) to include in 
#'   calculation of iqlr (smaller bounds means more stringent inclusion criteria)
#'
#' @description Takes idea from Wu et al. (citation below) and calculates IQLR for 
#'   Lambda, potentially useful if you believe there is an invariant group of 
#'   categories (e.g., taxa / genes) that are not changing (in absolute abundance) 
#'   between samples. IQLR is defined as 
#'   \deqn{IQLR_x = log(x_i/g(IQVF))}
#'   for i in 1,...,D. 
#'   IQVF are the CLR coordinates whose variance is within the inter-quantile range
#'   (defined by \code{probs} argument to this function). 
#'   A different IQVF is fit for each posteior sample as the IQVFs are calculted
#'   based on posterior estimates for Lambda. The variance of a CLR coordinate
#'   is defined as the norm of each row of Lambda[,focus.cov] (i.e., 
#'   the covariation in Eta, explained by those covariates). This definition of 
#'   variance allows uses to exclude variation from technical / trivial sources
#'   in calculation of IQVF/IQLR.   
#' 
#' @export
#' @return array of dimension (D, Q, iter) where D is number of taxa, Q is number
#' of covariates, and iter is number of posterior samples. 
#' 
#' @examples 
#' sim <- pibble_sim()
#' fit <- pibble(sim$Y, sim$X)
#' # Use first two covariates to define iqlr, just show first 5 samples
#' lambda_to_iqlr(fit, 1:2)[,,1:5] 
#' 
#' @references Jia R. Wu, Jean M. Macklaim, Briana L. Genge, Gregory B. Gloor (2017)
#'   Finding the center: corrections for asymmetry in high-throughput sequencing
#'   datasets. arxiv:1704.01841v1
lambda_to_iqlr <- function(m, focus.cov=NULL, probs=c(.25, .75)){
  req(m, "Lambda") # defensive
  if (!is.null(focus.cov)) focus.cov <- 1:m$Q
  if (is.character(focus.cov)) focus.cov <- which(focus.cov %in% names_categories)
  in.iqr <- matrix(0, ncategories(m), niter(m))
  
  # Convert to clr for calculating inter-quartile variable features (iqvf)
  m <- to_clr(m)
  
  # Calculate quantiles based on vector magnitude / covariance
  for (i in 1:niter(m)){
    L <- vec_to_mat(m$Lambda[,focus.cov,i])
    L <- rowSums(L*L)
    q <- stats::quantile(L, probs)
    in.iqr[,i] <- (L>= q[1]) & (L <= q[2]) 
  }
  
  # Transform 
  m <- to_proportions(m)
  Lambda_iqlr <- array(0, dim=c(ncategories(m), ncovariates(m), niter(m)))
  for (i in 1:niter(m)){
    V <- matrix(0, ncategories(m), ncategories(m))
    iqvf <- in.iqr[,i] ==1
    V[,iqvf] <- 1/sum(iqvf)
    V <- diag(nrow(V)) - V
    Lambda_iqlr[,,i] <- V%*%log(m$Lambda[,,i])
  }
  return(Lambda_iqlr)
}
