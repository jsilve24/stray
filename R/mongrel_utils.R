#' Provide random initialization for mongrel model 
#' 
#' Randomly initializes based on ALR transform of counts
#' plus random pseudocounts uniformily distributed between 
#' 0 and 1. 
#' 
#' Notation: \code{N} is number of samples and
#' \code{D} is number of multinomial categories
#'
#' @param Y matrix (D x N) of counts
#'
#' @return (D-1) x N matrix 
#' @importFrom driver alr
#' @export
#'
#' @examples
#' Y <- matrix(sample(1:100, 100), 10, 10)
#' random_mongrel_init(Y)
random_mongrel_init <- function(Y){
  N <- ncol(Y)
  D <- nrow(Y)
  t(driver::alr(t(Y)+runif(N*D)))
}


#' Check vector/matrix/data.frame for expected dimensions or throw error
#'
#' @param x object to check
#' @param d expected dimensions
#' @param par character name of x (for error message)
#'
#' @return nothing if no error, otherwise throws error
#' @export
#'
#' @examples
#' y <- c(1,3,4)
#' check_dims(y, 3, "y")
check_dims <- function(x, d, par){
  if (is.vector(x)){
    dc <- length(x)
  } else {
    dc <- dim(x)
  }
  d <- as.integer(d)
  
  if (!identical(dc, d)){
    s1 <- paste(d, collapse = " x ")
    s2 <- paste(dc, collapse = " x")
    s <- paste("Parameter", par, "should have dimension", s1, 
               "instead found dimensions", s2)
    stop(s)
  }
}

