#' Create pibblefit object
#'
#' @param D number of multinomial categories 
#' @param N number of samples
#' @param Q number of covariates
#' @param iter number of posterior samples
#' @param coord_system coordinate system objects are represented in (options 
#'   include "alr", "clr", "ilr", and "proportions")
#' @param alr_base integer category used as reference 
#'   (required if coord_system=="alr")
#' @param ilr_base (D x D-1) contrast matrix (required if coord_system=="ilr")
#' @param Eta Array of samples of Eta
#' @param Lambda Array of samples of Lambda
#' @param Sigma Array of samples of Sigma (null if coord_system=="proportions")
#' @param Sigma_default Array of samples of Sigma in alr base D, used if 
#'   coord_system=="proportions"
#' @param Y DxN matrix of observed counts
#' @param X QxN design matrix
#' @param upsilon scalar prior dof of inverse wishart prior
#' @param Theta prior mean of Lambda
#' @param Xi Matrix of prior covariance for inverse wishart 
#'   (null if coord_system=="proportions")
#' @param Xi_default Matrix of prior covariance for inverse wishart in alr 
#'   base D (used if coord_system=="proportions")
#' @param Gamma QxQ covariance matrix prior for Lambda
#' @param init matrix initial guess for Lambda used for optimization
#' @param names_categories character vector
#' @param names_samples character vector
#' @param names_covariates character vector
#'
#' @return object of class pibblefit
#'
#' @export 
#' @seealso \code{\link{pibble}}
pibblefit <- function(D, N, Q, coord_system, iter=NULL,  
                       alr_base=NULL, ilr_base=NULL,
                       Eta=NULL, Lambda=NULL, Sigma=NULL, Sigma_default=NULL, 
                       Y=NULL, X=NULL, upsilon=NULL, 
                       Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                       init=NULL, names_categories=NULL, names_samples=NULL, 
                       names_covariates=NULL){
  m <- new_pibblefit(D, N, Q, coord_system, iter, alr_base, ilr_base,
                      Eta, Lambda, Sigma, Sigma_default, 
                      Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                      init, names_categories, names_samples, 
                      names_covariates)
  verify(m)
  return(m)
}


new_pibblefit <- function(D, N, Q, coord_system, iter=NULL, 
                           alr_base=NULL, ilr_base=NULL,
                           Eta=NULL, Lambda=NULL,Sigma=NULL, Sigma_default=NULL, 
                           Y=NULL, X=NULL, upsilon=NULL, 
                           Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                           init=NULL, names_categories=NULL, names_samples=NULL, 
                           names_covariates=NULL){

  structure(
    list(
      # Basic Dimensions
      D=D, N=N, Q=Q, iter=iter, 
      # coordinate system parameters
      coord_system=coord_system, alr_base=alr_base, ilr_base=ilr_base, 
      # Results
      Eta=Eta, Lambda=Lambda, Sigma=Sigma, Sigma_default=Sigma_default, 
      # Data
      Y=Y, X=X, 
      # Priors
      upsilon=upsilon, Theta=Theta, Xi=Xi, Xi_default=Xi_default, Gamma=Gamma, 
      # Other
      init=init, names_categories=names_categories, names_samples=names_samples, 
      names_covariates=names_covariates
    ),
    class=c("pibblefit")
  )
}



#' Create orthusfit object
#'
#' @param D number of multinomial categories 
#' @param N number of samples
#' @param Q number of covariates
#' @param P Dimension of second dataset (e.g., nrows(Z) )
#' @param iter number of posterior samples
#' @param coord_system coordinate system objects are represented in (options 
#'   include "alr", "clr", "ilr", and "proportions")
#' @param alr_base integer category used as reference 
#'   (required if coord_system=="alr")
#' @param ilr_base (D x D-1) contrast matrix (required if coord_system=="ilr")
#' @param Eta Array of samples of Eta
#' @param Lambda Array of samples of Lambda
#' @param Sigma Array of samples of Sigma (null if coord_system=="proportions")
#' @param Sigma_default Array of samples of Sigma in alr base D, used if 
#'   coord_system=="proportions"
#' @param Y DxN matrix of observed counts
#' @param X QxN design matrix
#' @param Z PxN matrix of real valued observations
#' @param upsilon scalar prior dof of inverse wishart prior
#' @param Theta prior mean of Lambda
#' @param Xi Matrix of prior covariance for inverse wishart 
#'   (null if coord_system=="proportions")
#' @param Xi_default Matrix of prior covariance for inverse wishart in alr 
#'   base D (used if coord_system=="proportions")
#' @param Gamma QxQ covariance matrix prior for Lambda
#' @param init matrix initial guess for Lambda used for optimization
#' @param names_categories character vector
#' @param names_samples character vector
#' @param names_covariates character vector
#' @param names_Zdimensions character vector
#'
#' @return object of class pibblefit
#'
#' @export 
#' @seealso \code{\link{pibble}}
orthusfit <- function(D, N, Q, P, coord_system, iter=NULL,  
                      alr_base=NULL, ilr_base=NULL,
                      Eta=NULL, Lambda=NULL, Sigma=NULL, Sigma_default=NULL, 
                      Z=NULL, 
                      Y=NULL, X=NULL, upsilon=NULL, 
                      Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                      init=NULL, names_categories=NULL, names_samples=NULL, 
                      names_Zdimensions=NULL,
                      names_covariates=NULL){
  m <- new_orthusfit(D, N, Q, P, coord_system, iter, alr_base, ilr_base,
                     Eta, Lambda, Sigma, Sigma_default, 
                     Z,
                     Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                     init, names_categories, names_samples, 
                     names_Zdimensions,
                     names_covariates)
  verify(m)
  return(m)
}


new_orthusfit <- function(D, N, Q, P, coord_system, iter=NULL, 
                          alr_base=NULL, ilr_base=NULL,
                          Eta=NULL, Lambda=NULL,Sigma=NULL, Sigma_default=NULL, 
                          Z = NULL, 
                          Y=NULL, X=NULL, upsilon=NULL, 
                          Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                          init=NULL, names_categories=NULL, names_samples=NULL,
                          names_Zdimensions=NULL,
                          names_covariates=NULL){
  
  structure(
    list(
      # Basic Dimensions
      D=D, N=N, Q=Q, P=P, iter=iter, 
      # coordinate system parameters
      coord_system=coord_system, alr_base=alr_base, ilr_base=ilr_base, 
      # Results
      Eta=Eta, Lambda=Lambda, Sigma=Sigma, Sigma_default=Sigma_default, 
      # Data
      Y=Y, X=X, Z=Z,
      # Priors
      upsilon=upsilon, Theta=Theta, Xi=Xi, Xi_default=Xi_default, Gamma=Gamma, 
      # Other
      init=init, names_categories=names_categories, names_samples=names_samples, 
      names_Zdimensions=names_Zdimensions,
      names_covariates=names_covariates
    ),
    class=c("orthusfit")
  )
}

#' convert list to pibblefit
#' @param object list object
as.pibblefit <- function(object){
  class(object) <- "pibblefit"
  verify(object)
  return(object)
}

#' convert list to orthusfit
#' @param object list object
as.orthusfit <- function(object){
  class(object) <- "orthusfit"
  verify(object)
  return(object)
}


# simple function
# if x is NOT null evaluate expression y otherwise do nothing
ifnotnull <- function(x, y){
  if(!is.null(x)){
    return(eval(quote(y)))
  }
  else(return(x))
}

#' Simple verification of passed pibblefit object
#' @param m an object of class pibblefit
#' @param ... not used
#' @return throws error if any verification tests fail
#' @export 
verify.pibblefit <- function(m,...){
  # check basic dimensions that must always be present
  stopifnot(is.integer(m$N), is.integer(m$Q),
            is.integer(m$D))
  stopifnot(is.character(m$coord_system))
  if (m$coord_system == "ilr") stopifnot(!is.null(m$ilr_base))
  if (m$coord_system == "alr") stopifnot(!is.null(m$alr_base))
  if (m$coord_system != "proportions") stopifnot(is.null(m$Sigma_default))
  if (m$coord_system != "proportions") stopifnot(is.null(m$Xi_default))
  
  N <- m$N; D <- m$D; Q <- m$Q; iter <- m$iter
  Dm1 <- ifelse (m$coord_system %in% c("ilr", "alr"), D-1, D)
  
  # throw error if iter is null but Eta, Sigma, and Lambda are not. 
  if (is.null(m$iter)) stopifnot(is.null(m$Eta), is.null(m$Lambda), 
                                 is.null(m$Sigma), is.null(m$Sigma_default))
  ifnotnull(m$iter, stopifnot(is.integer(m$iter)))
  ifnotnull(m$Eta,check_dims(m$Eta, c(Dm1, N, iter), "pibblefit param Eta"))
  ifnotnull(m$Lambda,check_dims(m$Lambda, c(Dm1, Q, iter), "pibblefit param Lambda"))
  ifnotnull(m$Sigma,check_dims(m$Sigma, c(Dm1, Dm1, iter), "pibblefit param Sigma"))
  ifnotnull(m$Sigma_default,check_dims(m$Sigma_default, c(D-1, D-1, iter), "pibblefit param Sigma_default"))
  ifnotnull(m$Y,check_dims(m$Y, c(D, N), "pibblefit param Y"))
  ifnotnull(m$X,check_dims(m$X, c(Q, N), "pibblefit param X"))
  ifnotnull(m$upsilon,check_dims(m$upsilon, c(1), "pibblefit param upsilon"))
  ifnotnull(m$Theta,check_dims(m$Theta, c(Dm1, Q), "pibblefit param Theta"))
  ifnotnull(m$Xi,check_dims(m$Xi, c(Dm1,Dm1), "pibblefit param Xi"))
  ifnotnull(m$Xi_default,check_dims(m$Xi_default, c(D-1, D-1), "pibblefit param Xi_default"))
  ifnotnull(m$Gamma,check_dims(m$Gamma, c(Q, Q), "pibblefit param Gamma"))
  ifnotnull(m$init,check_dims(m$init, c(Dm1, N), "pibblefit param init"))
  ifnotnull(m$names_categories,check_dims(m$names_categories, c(D), "pibblefit param names_categories"))
  ifnotnull(m$names_samples,check_dims(m$names_samples, c(N), "pibblefit param names_samples"))
  ifnotnull(m$names_covariates,check_dims(m$names_covariates, c(Q), "pibblefit param names_covariates"))
}

#' Simple verification of passed orthusfit object
#' @param m an object of class orthusfit
#' @param ... not used
#' @return throws error if any verification tests fail
#' @export 
verify.orthusfit <- function(m,...){
  # check basic dimensions that must always be present
  stopifnot(is.integer(m$N), is.integer(m$Q),
            is.integer(m$D), is.integer(m$P))
  stopifnot(is.character(m$coord_system))
  if (m$coord_system == "ilr") stopifnot(!is.null(m$ilr_base))
  if (m$coord_system == "alr") stopifnot(!is.null(m$alr_base))
  if (m$coord_system != "proportions") stopifnot(is.null(m$Sigma_default))
  if (m$coord_system != "proportions") stopifnot(is.null(m$Xi_default))
  
  N <- m$N; D <- m$D; Q <- m$Q; iter <- m$iter; P <- m$P
  Dm1 <- ifelse (m$coord_system %in% c("ilr", "alr"), D-1, D)
  
  # throw error if iter is null but Eta, Sigma, and Lambda are not. 
  if (is.null(m$iter)) stopifnot(is.null(m$Eta), is.null(m$Lambda), 
                                 is.null(m$Sigma), is.null(m$Sigma_default))
  ifnotnull(m$iter, stopifnot(is.integer(m$iter)))
  ifnotnull(m$Eta,check_dims(m$Eta, c(Dm1, N, iter), "orthusfit param Eta"))
  ifnotnull(m$Lambda,check_dims(m$Lambda, c(Dm1+P, Q, iter), "orthusfit param Lambda"))
  ifnotnull(m$Sigma,check_dims(m$Sigma, c(Dm1+P, Dm1+P, iter), "orthusfit param Sigma"))
  ifnotnull(m$Sigma_default,check_dims(m$Sigma_default, c(D-1+P, D-1+P, iter), "orthusfit param Sigma_default"))
  ifnotnull(m$Y,check_dims(m$Y, c(D, N), "orthusfit param Y"))
  ifnotnull(m$Z,check_dims(m$Z, c(P, N), "orthusfit param Z"))
  ifnotnull(m$X,check_dims(m$X, c(Q, N), "orthusfit param X"))
  ifnotnull(m$upsilon,check_dims(m$upsilon, c(1), "orthusfit param upsilon"))
  ifnotnull(m$Theta,check_dims(m$Theta, c(Dm1+P, Q), "orthusfit param Theta"))
  ifnotnull(m$Xi,check_dims(m$Xi, c(Dm1+P,Dm1+P), "orthusfit param Xi"))
  ifnotnull(m$Xi_default,check_dims(m$Xi_default, c(D-1+P, D-1+P), "orthusfit param Xi_default"))
  ifnotnull(m$Gamma,check_dims(m$Gamma, c(Q, Q), "orthusfit param Gamma"))
  ifnotnull(m$init,check_dims(m$init, c(Dm1, N), "orthusfit param init"))
  ifnotnull(m$names_categories,check_dims(m$names_categories, c(D), "orthusfit param names_categories"))
  ifnotnull(m$names_Zdimensions,check_dims(m$names_Zdimensions, c(P), "orthusfit param names_Zdimensions"))
  ifnotnull(m$names_samples,check_dims(m$names_samples, c(N), "orthusfit param names_samples"))
  ifnotnull(m$names_covariates,check_dims(m$names_covariates, c(Q), "orthusfit param names_covariates"))
}

#' require elements to be non-null in pibblefit or throw error
#' @inheritParams req
#' @export 
req.pibblefit<- function(m, r){
  present <- sapply(m[r], is.null)
  if(any(present)){
    stop("pibblefit object does not contain required components:", r[present])
  }
}

#' require elements to be non-null in orthusfit or throw error
#' @inheritParams req
#' @export 
req.orthusfit<- function(m, r){
  req.pibblefit(m,r) # this works here
}

