#' Create maltipoofit object 
#' 
#' @inheritParams pibblefit
#' @inheritParams maltipoo_fit
#' @param VCScale scale factors (delta) for variance components 
#' @return object of class maltipoofit
#' @export
#' @seealso \code{\link{maltipoo}}
maltipoofit <- function(D, N, Q, P, coord_system, iter=NULL,  
                        alr_base=NULL, ilr_base=NULL,
                        Eta=NULL, Lambda=NULL, Sigma=NULL, Sigma_default=NULL, 
                        Y=NULL, X=NULL, upsilon=NULL, 
                        Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                        init=NULL, deltainit=NULL, names_categories=NULL, names_samples=NULL, 
                        names_covariates=NULL, VCScale=NULL, U=NULL){
  m <- new_maltipoo(D, N, Q, coord_system, iter, alr_base, ilr_base,
                    Eta, Lambda, Sigma, Sigma_default, 
                    Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                    init, deltainit, names_categories, names_samples, 
                    names_covariates, VCScale, U)
  verify(m)
  return(m)
}

# internal function 
new_maltipoofit <- function(D, N, Q, P, coord_system, iter=NULL, 
                            alr_base=NULL, ilr_base=NULL,
                            Eta=NULL, Lambda=NULL,Sigma=NULL, Sigma_default=NULL, 
                            Y=NULL, X=NULL, upsilon=NULL, 
                            Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                            init=NULL, deltainit=NULL, names_categories=NULL, names_samples=NULL, 
                            names_covariates=NULL, VCScale=NULL, U=NULL){
  m <- new_pibblefit(D, N, Q, coord_system, iter, alr_base, ilr_base,
                      Eta, Lambda, Sigma, Sigma_default, 
                      Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                      init, deltainit, names_categories, names_samples, 
                      names_covariates)
  m$VCScale <- VCScale
  m$U <- U
  m$deltainit <- deltainit
  m$P <- P
  class(m) <- c("maltipoofit", "pibblefit")
}

#' Simple verification of passed multipoo object
#' @param m an object of class multipoo
#' @param ... not used
#' @return throws error if any verification tests fail
#' @export 
verify.multipoo <- function(m,...){
  verify.pibblefit(m)
  stopifnot(is.integer(m$P))
  P <- m$P
  ifnotnull(m$VCScale, check_dims(m$VCScale, P))
  ifnotnull(m$U, check_dims(m$U, P*Q, Q))
  ifnotnull(m$deltainit, chec_dims(m$P))
}

#' require elements to be non-null in pibblefit or throw error
#' @inheritParams req
#' @export 
req.maltipoo<- function(m, r){
  present <- sapply(m[r], is.null)
  if(any(present)){
    stop("maltipoofit object does not contain required components:", r[present])
  }
}