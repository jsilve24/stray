#' Provide random initialization for pibble model 
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
#' random_pibble_init(Y)
random_pibble_init <- function(Y){
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
    dc <- unname(dim(x))
  }
  d <- unname(as.integer(d))
  
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
  if (is.null(m$names_categories)) return(NULL)
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

# For orthus objects
assign_combo_names <- function(m){
  # name-foo for some names of ortho
  if (!is.null(m$names_categories)){
    cnames <- assign_cat_names(m)
  } else {
    cnames <- paste0("c", 1:ncategories(m))
  }
  if (!is.null(m$names_Zdimensions)){
    znames <- m$names_Zdimensions
  } else {
    znames <- paste0("z", 1:ncategories(m))
  }
  combonames <- c(cnames, znames)
  return(combonames)
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
    } else if (identical(dimvars[[i]], "combo")) {
      n[[i]] <- assign_combo_names(m)
    } else if (is.data.frame(X)){
      if (max(X[,names(dimvars)[i]], na.rm=TRUE) == length(dimvars[[i]])) {
        n[[i]] <- dimvars[[i]]
      } else {
        n[i] <- list(NULL)
      }
    } else if (length(dimvars[[i]]) == dim(X)[i]){
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
name_array <- function(X, m, dimvars){
  n <- apply_names(X, m, dimvars)
  dimnames(X) <- n
  return(X)
}

# same as apply_names array but dimvars list must be 
# named with colnames of X to replace/use 
# @param as_factor if TRUE returns names as factors
name_tidy <- function(X, m, dimvars, as_factor=FALSE){
  if (is.null(names(dimvars))) stop("list element of dimvars must be named")
  n <- apply_names(X, m, dimvars)
  for (i in seq_along(dimvars)){
    d <- names(dimvars)[i]
    if (as_factor){
      if (!is.null(n[[i]])) X[[d]] <- factor(n[[i]][X[[d]]], levels=n[[i]])
    } else {
      if (!is.null(n[[i]])) X[[d]] <- n[[i]][X[[d]]]  
    }
  }
  return(X)
}

#' S3 for pibblefit apply names to pibblefit object
#' @param m object of class pibblefit
#' @param ... currently ignored
#' @return object of class pibblefit
name.pibblefit <- function(m, ...){
  if (!is.null(m$Eta)) {
    m$Eta <- name_array(m$Eta, m, list("cat", "sam", NULL))
  }
  if (!is.null(m$Lambda)){
    m$Lambda <- name_array(m$Lambda, m, list("cat", "cov", NULL))    
  }
  if (!is.null(m$Sigma)){
    m$Sigma <- name_array(m$Sigma, m, list("cat", "cat", NULL))
  } else if (!is.null(m$Sigma_default)) {
    m$Sigma_default <- name_array(m$Sigma_default, m,
                                         list("cat", "cat", NULL))
  }
  if (!is.null(m$Xi)){
    m$Xi <- name_array(m$Xi, m, list("cat", "cat"))
  } else if (!is.null(m$Xi_default)) { 
    m$Xi_default <- name_array(m$Xi_default, m, list("cat", "cat"))
  }
  if (!is.null(m$Theta)) {
    m$Theta <- name_array(m$Theta, m, list("cat", "cov"))
  }
  if (!is.null(m$Gamma)){
    m$Gamma <- name_array(m$Gamma, m, list("cov", "cov"))
    
  }
  if (!is.null(m$init)){
    m$init <- name_array(m$init, m, list("cat", "sam"))
  }
  return(m)
}


#' S3 for orthusfit apply names to orthusfit object
#' 
#' To avoid confusion, assigned default names to multinomial categories (c1 etc...) 
#' and zdimensions (z1 etc...)
#' 
#' @param m object of class orthusfit
#' @param ... currently ignored
#' @return object of class orthusfit
#' @export
name.orthusfit <- function(m, ...){
  if (!is.null(m$Eta)) {
    m$Eta <- name_array(m$Eta, m, list("cat", "sam", NULL))
  }

  combonames <- assign_combo_names(m)
  
  if (!is.null(m$Lambda)){
    m$Lambda <- name_array(m$Lambda, m, list(NULL, "cov", NULL))    
    dimnames(m$Lambda)[[1]] <- combonames
  }
  l <- list(combonames, combonames, NULL)
  if (!is.null(m$Sigma)){
    dimnames(m$Sigma) <- l
  } else if (!is.null(m$Sigma_default)) {
    dimnames(m$Sigma_default) <- l
  }
  if (!is.null(m$Xi)){
    dimnames(m$Xi) <- l[1:2]
  } else if (!is.null(m$Xi_default)) { 
    dimnames(m$Xi_default) <- l[1:2]
  }
  if (!is.null(m$Theta)) {
    m$Theta <- name_array(m$Theta, m, list(NULL, "cov"))
    dimnames(m$Theta)[[1]] <- combonames
  }
  if (!is.null(m$Gamma)){
    m$Gamma <- name_array(m$Gamma, m, list("cov", "cov"))
  }
  if (!is.null(m$init)){
    m$init <- name_array(m$init, m, list("cat", "sam"))
  }
  return(m)
}


# ifelse that can handle vector/matrix/dataframe/array output - inline but not 
# vectorized
mifelse <- function(b, y, n){
  if (b) return(y) else return(n)
}


# x is vector of elements that should either be NULL or all the same 
try_set_dims <- function(x){
  if (length(x) ==0){
    stop("Not Information provided to set dimension")
  }
  if (all(x[1]==x)) {
    return(as.integer(x[1])) 
  }else{
    msg <- paste(x, collapse = ",")
    msg <- paste("Dimension missmatch in arguments: [", msg, "]", sep="")
    stop(msg)
  } 
}


# to make default checking of extra arguments easier
args_null <- function(par, argl, default){
  if (is.null(argl[[par]])) return(default)
  return(argl[[par]])
}

# Parse Timer
parse_timer_seconds <- function(timer){
  if (is.null(timer)) return(NULL)
  n <- strsplit(names(timer), "_")
  n1 <- unlist(lapply(n, FUN = function(x) x[1]))
  n1n <- as.integer(as.factor(n1))
  n1s <- split(timer, n1n)
  n1names <- lapply(split(n1, n1n), unique)
  times <- lapply(n1s, diff)
  times.seconds <- unlist(times)/1e9
  names(times.seconds) <- unlist(n1names)
  return(times.seconds)
}



