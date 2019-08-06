
#' Interface to fit basset models
#' 
#' Basset (A Lazy Learner) - non-linear regression models in stray
#'
#' @param Y D x N matrix of counts (if NULL uses priors only)
#' @param X Q x N matrix of covariates (cannot be NULL)
#' @param upsilon dof for inverse wishart prior (numeric must be > D) 
#'   (default: D+3)
#' @param Theta A function from dimensions dim(X) -> (D-1)xN (prior mean of gaussian process)
#' @param Gamma A function from dimension dim(X) -> NxN (kernel matrix of gaussian process)
#' @param Xi (D-1)x(D-1) prior covariance matrix
#'   (default: ALR transform of diag(1)*(upsilon-D)/2 - this is 
#'   essentially iid on "base scale" using Aitchison terminology)
#' @param init (D-1) x Q initialization for Eta for optimization 
#' @param pars character vector of posterior parameters to return
#' @param m object of class bassetfit 
#' @param ... other arguments passed to \link{pibble} (which is used internally to 
#'  fit the basset model)
#' 
#' @details the full model is given by:
#'    \deqn{Y_j \sim Multinomial(Pi_j)} 
#'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
#'    \deqn{Eta \sim MN_{D-1 x N}(Lambda, Sigma, I_N)}
#'    \deqn{Lambda \sim GP_{D-1 x Q}(Theta(X), Sigma, Gamma(X))}
#'    \deqn{Sigma \sim InvWish(upsilon, Xi)}
#'  Where Gamma(X) is short hand for the Gram matrix of the Kernel function. 
#'  
#'  Default behavior is to use MAP estimate for uncollaping the LTP 
#'  model if laplace approximation is not preformed. 
#' @return an object of class bassetfit
#' @md
#' @name basset_fit
NULL

#' @rdname basset_fit
#' @export
basset <- function(Y=NULL, X, upsilon=NULL, Theta=NULL, Gamma=NULL, Xi=NULL, 
                   init=NULL, pars=c("Eta", "Lambda", "Sigma"), ...){
  
  if (!is.null(Theta)) {
    Theta_train <- Theta(X)
  } else {
    Theta <- function(X) matrix(0, nrow(Y)-1, ncol(X)) 
    Theta_train <- Theta(X)
  }
  if (!is.null(Gamma)) {
    Gamma_train <- Gamma(X)
  } else {
    stop("No Default Kernel For Gamma Implemented")
  }
  
  out <- pibble(Y, X=diag(ncol(X)), upsilon, Theta_train, Gamma_train, Xi, init, pars, ...)
  out$Q <- as.integer(nrow(X))
  out$X <- X
  out$Theta <- Theta
  out$Gamma <- Gamma
  class(out) <- c("bassetfit", "pibblefit")
  verify(out)
  return(out)
}

#' @rdname basset_fit
#' @export
refit.bassetfit <- function(m, pars=c("Eta", "Lambda", "Sigma"), ...){
  # Store coordinates and tranfsorm to cannonical representation
  l <- store_coord(m)
  m <- to_alr(m, m$D)
  
  # Concatenate parameters to pass to basset function
  argl <- list(...)
  argl$pars <- pars
  ml <- as.list(m)
  argl <- c(ml, argl)
  
  # Need to handle iter as part of m but no n_samples passed
  # in this situation should pull iter from m and pass as n_samples to pibble 
  if (is.null(argl[["n_samples"]]) & !is.null(m$iter)) argl[["n_samples"]] <- m$iter 
  
  # pass to basset function
  m <- do.call(basset, argl)
  
  # Reapply original coordinates
  m <- reapply_coord(m, l)
  verify(m)
  return(m)
}