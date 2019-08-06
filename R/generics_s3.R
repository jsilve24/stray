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
sample_prior <- function(m, n_samples=2000L, ...){
  UseMethod("sample_prior", m)
}

#' Generic method for fitting model from passed model fit object
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
#' @param m An object of class pibblefit 
#' @details An alternative approach to accessing these dimensions is to 
#'   access them directly from the pibblefit object using list indexing. 
#' * \code{ncategories} is equivalent to \code{m$D}
#' * \code{nsamples} is equivalent to \code{m$N}
#' * \code{ncovariates} is equivalent to \code{m$Q}
#' @return integer 
#' @name access_dims
#' @examples 
#' \dontrun{
#' m <- pibble(Y, X)
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

#' @rdname access_dims
#' @export
niter <- function(m){
  UseMethod("niter", m)
}


#' Generic method for getting and setting dimension names of fit object
#' 
#' @param m object
#' @param value character vector (or NULL)
#' @name name_dims
#' 
#' @details \code{names_coords} is different than \code{names_categories}. 
#' \code{names_categories} provides access to the basic names of each multinomial 
#' category. In contrast, \code{names_coords} provides access to the 
#' names of the coordinates in which an object is represented. These coordinate
#' names are based on the category names. For example, category names may be, 
#' (OTU1, ..., OTUD) where as coordinate names could be (log(OTU1/OTUD), etc...)
#' if object is in default coordinate system. 
NULL

#' @rdname name_dims
#' @export
names_covariates <- function(m){
  UseMethod("names_covariates", m)
}

#' @rdname name_dims
#' @export
names_samples <- function(m){
  UseMethod("names_samples", m)
}


#' @rdname name_dims
#' @export
names_categories <- function(m){
  UseMethod("names_categories", m)
}

#' @rdname name_dims
#' @export
names_coords <- function(m){
  UseMethod("names_coords", m)
}


#' @rdname name_dims
#' @export
`names_covariates<-` <- function(m, value){
  UseMethod("names_covariates<-", m)
}

#' @rdname name_dims
#' @export
`names_samples<-` <- function(m, value){
  UseMethod("names_samples<-", m)
}


#' @rdname name_dims
#' @export
`names_categories<-` <- function(m, value){
  UseMethod("names_categories<-", m)
}



#' Generic Method to Plot Posterior Predictive Summaries 
#' 
#' @param m model object
#' @param ... other arguments to pass 
#' 
#' @return vector
#' @name ppc_summary
NULL

#' @rdname ppc_summary
#' @export
ppc_summary <- function(m, ...){
  UseMethod("ppc_summary", m)
}



