#' Transform Fit Mongrel Parameters to other representations
#' 
#' These are a collection of convenience functions for transforming
#' mongrel fit objects to a number of different representations includeing
#' ILR bases, CLR coordinates, ALR coordinates, and proportions. 
#'
#' @param m object of class mongrelfit (e.g., output of \code{\link{mongrel}})
#' @param d (integer) multinomial category to take as new alr reference
#' @param V (matrix) contrast matrix for ILR basis to transform into to
#'
#' @details Note that there is a degeneracy of representations for a covariance 
#' matrix represented in terms of proportions. As such the function 
#' \code{mongrel_to_proportions} does not attempt to transform parameters Sigma
#' or prior Xi and instead just remove them from the mongrelfit object returned
#'
#' @return mongrelfit object
#' @export
#' @name mongrel_transforms
NULL

#' @rdname mongrel_transforms
#' @export
mongrel_to_alr <- function(m, d){
  if (m$coord_system != "default") {
    stop("should be applied to original output of mongrel, this object has already been transformed")
  }
  m$Eta <- alrInv_array(m$Eta, m$D, 1)
  m$Eta <- alr_array(m$Eta, d, 1)
  m$Lambda <- alrInv_array(m$Lambda, m$D, 1)
  m$Lambda <- alr_array(m$Lambda, d, 1)
  for (i in 1:m$iter){
    m$Sigma[,,i] <- alrvar2alrvar(m$Sigma[,,i], m$D, d)
  }
  # Transform priors as well 
  m$Xi <- alrvar2alrvar(m$Xi, m$D, d)
  m$Theta <- alrInv_array(m$Theta, m$D, 1)
  m$Theta <- alr_array(m$Theta, d, 1)
  
  m$summary <- NULL
  m$coord_system <- "alr"
  return(m)
}

#' @rdname mongrel_transforms
#' @export
mongrel_to_ilr <- function(m, V){
  if (m$coord_system != "default") {
    stop("should be applied to original output of mongrel, this object has already been transformed")
  }
  m$Eta <- alrInv_array(m$Eta, m$D, 1)
  m$Eta <- ilr_array(m$Eta, V, 1)
  m$Lambda <- alrInv_array(m$Lambda, m$D, 1)
  m$Lambda <- ilr_array(m$Lambda, V, 1)
  for (i in 1:m$iter){
    m$Sigma[,,i] <- alrvar2ilrvar(m$Sigma[,,i], m$D, V)
  }
  # Transform priors as well 
  m$Xi <- alrvar2ilrvar(m$Xi, m$D, V)
  m$Theta <- alrInv_array(m$Theta, m$D, 1)
  m$Theta <- ilr_array(m$Theta, V, 1)
  
  m$summary <- NULL
  m$coord_system <- "ilr"
  return(m)
}

#' @rdname mongrel_transforms
#' @export
mongrel_to_clr <- function(m){
  if (m$coord_system != "default") {
    stop("should be applied to original output of mongrel, this object has already been transformed")
  }
  m$Eta <- alrInv_array(m$Eta, m$D, 1)
  m$Eta <- clr_array(m$Eta, 1)
  m$Lambda <- alrInv_array(m$Lambda, m$D, 1)
  m$Lambda <- clr_array(m$Lambda, 1)
  Sigma <- array(0, dim=c(m$D, m$D, m$iter))
  for (i in 1:m$iter){
    Sigma[,,i] <- alrvar2clrvar(m$Sigma[,,i], m$D)
  }
  m$Sigma <- Sigma
  
  # Transform priors as well 
  m$Xi <- alrvar2clrvar(m$Xi, m$D)
  m$Theta <- alrInv_array(m$Theta, m$D, 1)
  m$Theta <- clr_array(m$Theta, 1)
  
  m$summary <- NULL
  m$coord_system <- "clr"
  return(m)
}

#' @rdname mongrel_transforms
#' @export
mongrel_to_proportions <- function(m){
  if (m$coord_system != "default") {
    stop("should be applied to original output of mongrel, this object has already been transformed")
  }
  m$Eta <- alrInv_array(m$Eta, m$D, 1)
  m$Lambda <- alrInv_array(m$Lambda, m$D, 1)
  m$Sigma <- NULL
  
  # Transform priors as well 
  m$Xi <- NULL
  m$Theta <- alrInv_array(m$Theta, m$D, 1)

  m$summary <- NULL
  m$coord_system <- "proportions"
  return(m)
} 