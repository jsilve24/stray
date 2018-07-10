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
  present <- sapply(m[r], is.null)
  if(any(present)){
    stop("object does not contain required components:", r[present])
  }
}