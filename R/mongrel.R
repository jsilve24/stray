#' mongrel: Fitting and Analysis of Multinomial Logistic Normal Regression Models
#' 
#'  Provides methods for fitting and inspection of Bayesian Multinomial 
#'  Logistic Normal Regression Models ("Mongrel Models") using MAP estimation 
#'  (with the ADAM optimizer) and Laplace Approximation. Key functionality is 
#'  implemented in C++ for scalability. 
#'  
#' @docType package
#' @name mongrel_package
#' 
#' @useDynLib mongrel
#' @importFrom Rcpp sourceCpp
NULL