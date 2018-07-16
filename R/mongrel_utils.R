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
  if (is.vector(x) | is.factor(x)){
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

# function to help with assignment of category naming conventions for different
# coordinate systems. 
assign_cat_names <- function(m){
  if (m$coord_system=="proportions") return(paste0("prop_", m$names_categories))
  else if (m$coord_system=="clr") return(paste0("clr_", m$names_categories))
  else if (m$coord_system=="alr"){
    n1 <- m$names_categories[-m$alr_base]
    n2 <- m$names_categories[m$alr_base]
    n <- paste0("log(",n1, "/", n2,")")
    return(n)
  } else if (m$coord_system=="ilr") return(colnames(m$ilr_base))
  else stop("not a recognized coordinate system to name")
}

apply_names <- function(X, m, dimvars){
  n <- list()
  for (i in seq_along(dimvars)){
    if (identical(dimvars[[i]], "sam") & !is.null(m$names_samples)){
      n[[i]] <- m$names_samples
    } else if (identical(dimvars[[i]], "cov") & !is.null(m$names_covariates)) {
      n[[i]] <- m$names_covariates
    } else if (identical(dimvars[[i]], "cat") & !is.null(m$names_categories)) {
      n[[i]] <- assign_cat_names(m)
    } else if (length(dimvars[[i]]) == dim(X)[i]) {
      n[[i]] <- dimvars[[i]]
    } else {
      n[i] <- list(NULL)
    }
  }
  if (!is.null(names(dimvars))) names(n) <- names(dimvars)
  return(n)
}

# dimvars = cat, sam, cov or NULL - list
# or can pass vector as element of the list 
apply_names_array <- function(X, m, dimvars){
  n <- apply_names(X, m, dimvars)
  dimnames(X) <- n
  return(X)
}

# same as apply_names array but dimvars list must be 
# named with colnames of X to replace/use 
apply_names_tidy <- function(X, m, dimvars){
  if (is.null(names(dimvars))) stop("list element of dimvars must be named")
  n <- apply_names(X, m, dimvars)
  for (i in seq_along(dimvars)){
    d <- names(dimvars)[i]
    if (!is.null(n[[i]])) X[[d]] <- n[[i]][X[[d]]]
  }
  return(X)
}

# apply names to mongrel object
apply_names_mongrel <- function(m){
  m$Eta <- apply_names_array(m$Eta, m, list("cat", "sam", NULL))
  m$Lambda <- apply_names_array(m$Lambda, m, list("cat", "cov", NULL))
  if (!is.null(m$Sigma)){
    m$Sigma <- apply_names_array(m$Sigma, m, list("cat", "cat", NULL))
    m$Xi <- apply_names_array(m$Xi, m, list("cat", "cat"))
  } else {
    m$Sigma_default <- apply_names_array(m$Sigma_default, m,
                                         list("cat", "cat", NULL))
    m$Xi_default <- apply_names_array(m$Xi_default, m, list("cat", "cat"))
  }
  m$Theta <- apply_names_array(m$Theta, m, list("cat", "cov"))
  m$Gamma <- apply_names_array(m$Gamma, m, list("cov", "cov"))
  m$init <- apply_names_array(m$init, m, list("cat", "sam"))
  return(m)
}