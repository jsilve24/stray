#' Multivariate RBF Kernel
#' 
#' Designed to be partially specified. (see examples)
#' @param X covariate (dimension Q x N; i.e., covariates x samples) 
#' @param sigma scalar parameter 
#' @param rho scalar bandwidth parameter
#' @param jitter small scalar to add to off-diagonal of gram matrix 
#'   (for numerical underflow issues)
#' @param c vector parameter defining intercept for linear kernel
#' 
#' @details Gram matrix G is given by 
#' 
#' SE (squared exponential):
#' \deqn{G = \sigma^2 * exp(-[(X-c)'(X-c)]/(s*\rho^2))}
#' 
#' LINEAR:
#' \deqn{G = \sigma^2*(X-c)'(X-c)}
#' 
#' 
#' @return Gram Matrix (N x N) (e.g., the Kernel evaluated at 
#' each pair of points)
#' @name kernels
#' @examples
#'   # Create Partial for use with basset
#'   K <- function(X) SE(X, 2, .2)
#'   
#'   # Example use
#'   X <- matrix(rnorm(10), 2, 5)
#'   G <- K(X)
#'   G # this is the gram matrix (the kernel evaluated on a finite set of points)
NULL

#' @rdname kernels
#' @export 
SE <- function(X, sigma=1, rho=median(as.matrix(dist(t(X)))), jitter=1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-dist^2/(2*rho^2)) + jitter*diag(ncol(dist))
  return(G)
}


#' @rdname kernels
#' @export
LINEAR <- function(X, sigma=1, c=rep(0, nrow(X))){
  E <- sweep(X, 1, c)
  G <- sigma^2*crossprod(E)
  return(G)
}