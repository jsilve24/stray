#' Generic method for verifying new objects
#' 
#' Intended to be called internally by package or object creator
#' 
#' @param m object
#' @param ... other arguments to be passed to verify
#' 
#' @return throws error if verify test fails
#' @export 
verify <- function(m, ...){
  UseMethod("verify", m)
}

#' Generic method for ensuring object contains required elements
#' 
#' Intended to be called internally by package 
#' 
#' @param m object
#' @param r vector of elements to test for 
#' 
#' @return throws error if required element is not present 
#' @export 
req <- function(m, r){
  UseMethod("req", m)
}


#' Generic method for applying names to an object
#' 
#' Intended to be called internally by package
#' 
#' @param m object
#' @param ... other arguments to be passed 
#' 
#' @return object of same class but with names applied to dimensions 
#' @export
name <- function(m, ...){
  UseMethod("name", m)
}


#' Generic method for sampling from prior distribution of object
#' 
#' @param m object
#' @param n_samples number of samples to produce
#' @param ... other arguments to be passed
#' 
#' @export
#' @return object of the same class 
sample_prior <- function(m, n_sample=2000, ...){
  UseMethod("sample_prior", m)
}

#' Generic method for fitting model from passed model fit objet
#' 
#' @param m object
#' @param ... other arguments passed that control fitting
#' 
#' @export
#' @return object of the same class as \code{m}
refit <- function(m, ...){
  UseMethod("refit", m)
}


#' Generic method for visualizing posterior predictive checks
#' @param m object
#' @param ... other arguments passed that control visualization
#' @export
ppc <- function(m, ...){
  UseMethod("ppc", m)
}


#' Generic method for accessing model fit dimensions
#'
#' @param m An object of class mongrelfit 
#' @details An alternative approach to accessing these dimensions is to 
#'   access them directly from the mongrelfit object using list indexing. 
#' * \code{ncategories} is equivalent to \code{m$D}
#' * \code{nsamples} is equivalent to \code{m$N}
#' * \code{ncovariates} is equivalent to \code{m$Q}
#' @return integer 
#' @name access_dims
#' @examples 
#' \dontrun{
#' m <- mongrel(Y, X)
#' ncategories(m)
#' nsamples(m)
#' ncovariates(m)
#' }
NULL

#' @rdname access_dims
#' @export
ncategories <- function(m){
  UseMethod("ncategories", m)
}

#' @rdname access_dims
#' @export
nsamples <- function(m){
  UseMethod("nsamples", m)
}

#' @rdname access_dims
#' @export
ncovariates <- function(m){
  UseMethod("ncovariates", m)
}

