#' Transform Fit Stray Parameters to other representations
#' 
#' These are a collection of convenience functions for transforming
#' stray fit objects to a number of different representations including
#' ILR bases, CLR coordinates, ALR coordinates, and proportions. 
#'
#' @param m object of class pibblefit (e.g., output of \code{\link{pibble}})
#' @param d (integer) multinomial category to take as new alr reference
#' @param V (matrix) contrast matrix for ILR basis to transform into to (defaults to 
#'   \code{driver::create_default_ilr_base(D)})
#'
#' @details Note: that there is a degeneracy of representations for a covariance 
#' matrix represented in terms of proportions. As such the function 
#' \code{to_proportions} does not attempt to transform parameters Sigma
#' or prior Xi and instead just removes them from the pibblefit object returned. 
#' 
#' @return pibblefit object
#' @name stray_transforms
#' @import driver 
#' @examples
#' \dontrun{
#' m <- pibble(Y, X)
#' m.prop <- to_proportions(m)
#' # convert back to default coordinates (alr with D-th part as reference)
#' m <- to_alr(m.prop, ncategories(m))
#' V <- driver::create_default_ilr_base(ncategories(m))
#' m.ilr <- to_ilr(m, V)
#' m.clr <- to_clr(m)
#' }
NULL


#' @rdname stray_transforms
#' @export
to_proportions <- function(m){
  if (m$coord_system == "alr"){
    if (!is.null(m$Eta)) m$Eta <- alrInv_array(m$Eta, m$alr_base, 1)
    if (!is.null(m$Lambda)) m$Lambda <- alrInv_array(m$Lambda, m$alr_base, 1)
    if (!is.null(m$Sigma)){
      if (m$alr_base != m$D){
        for (i in 1:m$iter){
          m$Sigma[,,i] <- alrvar2alrvar(m$Sigma[,,i], m$alr_base, m$D)
        }
      }
      m$Sigma_default <- m$Sigma
      m$Sigma <- NULL
    }
    # Transform Priors as well 
    if (!is.null(m$Xi)){
      if (m$alr_base != m$D){
        m$Xi <- alrvar2alrvar(m$Xi, m$alr_base, m$D)
      }
      m$Xi_default <- m$Xi
      m$Xi <- NULL
    }
    if (!is.null(m$Theta)){
      if (!inherits(m, "bassetfit")) m$Theta <- alrInv_array(m$Theta, m$alr_base, 1)
    }
    if (!is.null(m$init)) m$init <- alrInv_array(m$init, m$alr_base, 1)
  }
  if (m$coord_system == "ilr"){
    if (!is.null(m$Eta)) m$Eta <- ilrInv_array(m$Eta, m$ilr_base, 1)
    if (!is.null(m$Lambda)) m$Lambda <- ilrInv_array(m$Lambda, m$ilr_base, 1)
    if (!is.null(m$Sigma)){
      for (i in 1:m$iter){
        m$Sigma[,,i] <- ilrvar2alrvar(m$Sigma[,,i], m$ilr_base, m$D)
      }
      m$Sigma_default <- m$Sigma
      m$Sigma <- NULL
    }
    
    # Transform priors as well 
    if (!is.null(m$Xi)){
      m$Xi <- ilrvar2alrvar(m$Xi, m$ilr_base, m$D)
      m$Xi_default <- m$Xi
      m$Xi <- NULL  
    }
    if (!is.null(m$Theta)) {
      if (!inherits(m, "bassetfit")) m$Theta <- ilrInv_array(m$Theta, m$ilr_base, 1)  
    }
    if (!is.null(m$init)) m$init <- ilrInv_array(m$init, m$ilr_base, 1)
  }
  if (m$coord_system == "clr"){
    if (!is.null(m$Eta)) m$Eta <- clrInv_array(m$Eta, 1)
    if (!is.null(m$Lambda)) m$Lambda <- clrInv_array(m$Lambda, 1)
    if (!is.null(m$Sigma)){
      Sigma_default <- array(0, dim=c(m$D-1, m$D-1, m$iter))
      for (i in 1:m$iter){
        Sigma_default[,,i] <- clrvar2alrvar(m$Sigma[,,i], m$D)
      }
      m$Sigma <- NULL
      m$Sigma_default <- Sigma_default
    }
    # Transform priors as well
    if (!is.null(m$Xi)){
      m$Xi_default <- clrvar2alrvar(m$Xi, m$D)
      m$Xi <- NULL      
    }
    if (!is.null(m$Theta)){
      if (!inherits(m, "bassetfit")) m$Theta <- clrInv_array(m$Theta, 1)  
    }
    if (!is.null(m$init)) m$init <- clrInv_array(m$init, 1)
  }
  if (m$coord_system=="proportions"){
    return(m)
  }
  m$summary <- NULL
  m$coord_system <- "proportions"
  m$ilr_base <- NULL
  m$alr_base <- NULL
  return(m)
}




#' @rdname stray_transforms
#' @export
to_alr <- function(m, d){
  if (m$coord_system=="alr"){
    if (m$alr_base == d) return(m)
  }
  m <- to_proportions(m)
  
  if (!is.null(m$Eta)) m$Eta <- alr_array(m$Eta, d, 1)
  if (!is.null(m$Lambda)) m$Lambda <- alr_array(m$Lambda, d, 1)
  if (!is.null(m$Sigma)){
    m$Sigma <- array(0, dim=dim(m$Sigma_default))
    for (i in 1:m$iter){
      m$Sigma[,,i] <- alrvar2alrvar(m$Sigma_default[,,i], m$D, d)
    }
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- alrvar2alrvar(m$Xi_default, m$D, d)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)){
    if (!inherits(m, "bassetfit")) m$Theta <- alr_array(m$Theta, d, 1)  
  }
  if (!is.null(m$init)) m$init <- alr_array(m$init, d, 1)
  
  m$summary <- NULL
  m$coord_system <- "alr"
  m$alr_base <- d
  return(m)
}

#' @rdname stray_transforms
#' @export
to_ilr <- function(m, V=NULL){
  if (m$coord_system=="ilr"){
    if (all.equal(m$ilr_base, V)) return(m)
  }
  if (is.null(V)) V <- driver::create_default_ilr_base(m$D)
  m <- to_proportions(m)
  
  if (!is.null(m$Eta)) m$Eta <- ilr_array(m$Eta, V, 1)
  if (!is.null(m$Lambda)) m$Lambda <- ilr_array(m$Lambda, V, 1)
  if (!is.null(m$Sigma)){
    m$Sigma <- array(0, dim=dim(m$Sigma_default))
    for (i in 1:m$iter){
      m$Sigma[,,i] <- alrvar2ilrvar(m$Sigma_default[,,i], m$D, V)
    }
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- alrvar2ilrvar(m$Xi_default, m$D, V)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)){
    if (!inherits(m, "bassetfit")) m$Theta <- ilr_array(m$Theta, V, 1)  
  }
  if (!is.null(m$init)) m$init <- ilr_array(m$init, V, 1)
  
  m$summary <- NULL
  m$coord_system <- "ilr"
  m$ilr_base <- V
  return(m)
}

#' @rdname stray_transforms
#' @export
to_clr <- function(m){
  if (m$coord_system=="clr") return(m)
  m <- to_proportions(m)

  if (!is.null(m$Eta)) m$Eta <- clr_array(m$Eta, 1)
  if (!is.null(m$Lambda)) m$Lambda <- clr_array(m$Lambda, 1)
  if (!is.null(m$Sigma)){
    m$Sigma <- array(0, dim=c(m$D, m$D, m$iter))
    for (i in 1:m$iter){
      m$Sigma[,,i] <- alrvar2clrvar(m$Sigma_default[,,i], m$D)
    }
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- alrvar2clrvar(m$Xi_default, m$D)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)){
    if (!inherits(m, "bassetfit")) m$Theta <- clr_array(m$Theta, 1)  
  }
  if (!is.null(m$init)) m$init <- clr_array(m$init, 1)
  
  m$summary <- NULL
  m$coord_system <- "clr"
  return(m)
}

#' Holds information on coordinates system to later be reapplied
#' 
#' \code{store_coord} stores coordinate information for pibblefit object
#' and can be reapplied with function \code{reapply_coord}. Some coordinate
#' systems are not useful for computation and this makes it simple keep 
#' returned object from computations in the same coordinate system as the input. 
#' (Likely most useful inside of a package)
#' 
#' 
#' @param m object of class pibblefit
#' @param l object returned by function \code{store_coord}
#' @name store_coord
#' @return \code{store_coord} list with important information to identify c
#'  coordinate system of pibblefit object. \code{reapply_coord} pibblefit object
#'  in coordinate system previously stored. 
NULL


#' @rdname store_coord
#' @export
store_coord <- function(m){
  l <- list()
  l$coord_system <- m$coord_system
  l$alr_base <- m$alr_base
  l$ilr_base <- m$ilr_base
  return(l)
}

#' @rdname store_coord
#' @export
reapply_coord <- function(m, l){
  if (l$coord_system == "proportions") return(to_proportions(m))
  if (l$coord_system == "clr") return(to_clr(m))
  if (l$coord_system == "alr") return(to_alr(m, l$alr_base))
  if (l$coord_system == "ilr") return(to_ilr(m, l$ilr_base))
  stop("not a recognized coordinate system")
}