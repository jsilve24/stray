#' stray: Fitting and Analysis of Multinomial Logistic Normal  Models
#' 
#'  Provides methods for fitting and inspection of Bayesian Multinomial 
#'  Logistic Normal Models using MAP estimation 
#'  (with the ADAM optimizer) and Laplace Approximation. Key functionality is 
#'  implemented in C++ for scalability. 
#'  
#' @docType package
#' @name stray_package
#' 
#' @useDynLib stray
#' @importFrom Rcpp sourceCpp
NULL

globalVariables(".")



.onAttach <- function(libname, pkgname) {
  packageStartupMessage("The stray package has been renamed fido due to name collision with another package on CRAN.", "\n", "\n",
                        "At this time please switch to the fido package (devtools::install_github('jsilve24/fido')) where the project is now being actively developed.", "\n\n", 
                        "I appologize for the inconvenience - Justin Silverman")
}