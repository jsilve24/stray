#' Convert mongrel samples of Eta Lambda and Sigma to tidy format
#' 
#' Combines them all into a single tibble, see example for formating and 
#' column headers. Primarily designed to be used by 
#' \code{\link{summary.mongrelfit}}. 
#' 
#' @param m an object of class mongrelfit
#' @param use_names should dimension indicies be replaced by
#'   dimension names if provided in data used to fit mongrel model.  
#' 
#' @importFrom driver gather_array
#' @importFrom dplyr bind_rows group_by 
#' @export
#' @return tibble
#' @examples 
#' sim <- mongrel_sim()
#' fit <- mongrel(sim$Y, sim$X)
#' fit_tidy <- mongrel_tidy_samples(fit, use_names=TRUE)
#' head(fit_tidy)
mongrel_tidy_samples<- function(m, use_names=FALSE){
  l <- list()
  if (!is.null(m$Eta)) l$Eta <- driver::gather_array(m$Eta, val, coord, sample, iter)
  if (!is.null(m$Lambda)) l$Lambda <- driver::gather_array(m$Lambda, val, coord, covariate, iter)
  if (!is.null(m$Sigma)) l$Sigma <- driver::gather_array(m$Sigma, val, coord, coord2, iter) 

  l <- dplyr::bind_rows(l, .id="Parameter")
  
  if (!use_names) return(l)
  # Deal with names (first create conversion list)
  cl <- list()
  if (!is.null(m$Eta)) {
    cl[["sample"]] = "sam"
    cl[["coord"]] = "cat"
  }
  if (!is.null(m$Lambda)){
    cl[["covariate"]] = "cov"
    cl[["coord"]] = "cat"
  }
  if (!is.null(m$Sigma)){
    cl[["coord"]] = "cat"
    cl[["coord2"]] = "cat"
  }

  l <- name_tidy(l, m, cl)
  
  return(l)
}

# Internal function to check if summary has already been precomputed. 
summary_check_precomputed <- function(m, pars){
  if (!is.null(m$summary)){
    if (all(!is.null(m$summary[pars]))) return(TRUE)
  }
  return(FALSE)
}

#' Summarise mongrelfit object and print posterior quantiles
#' 
#' Default calculates median, mean, 50\% and 95\% credible interval
#' 
#' @param object an object of class mongrelfit 
#' @param pars character vector (default: c("Eta", "Lambda", "Sigma"))
#' @param use_names should summary replace dimension indicies with mongrelfit 
#'   names if names Y and X were named in call to \code{\link{mongrel}}
#' @param gather_prob if TRUE then prints quantiles in long format rather than 
#'  wide (useful for some plotting functions)
#' @param ... other expressions to pass to summarise (using name 'val' unquoted is 
#'   probably what you want)
#' @import dplyr
#' @importFrom driver summarise_posterior
#' @importFrom purrr map
#' @importFrom tidybayes mean_qi
#' @importFrom dplyr group_by select ungroup
#' @export
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' summary(fit, pars="Eta", median = median(val))
#' 
#' # Some later functions make use of precomputation
#' fit$summary <- summary(fit)
#' }
summary.mongrelfit <- function(object, pars=NULL, use_names=TRUE, gather_prob=FALSE,
                               ...){
  if (is.null(pars)) {
    pars <- c("Eta", "Lambda", "Sigma")
    pars <- pars[pars %in% names(object)] # only for the ones that are present 
  }
  
  
  # if already calculated
  if (summary_check_precomputed(object, pars)) return(object$summary[pars])
  
  mtidy <- dplyr::filter(mongrel_tidy_samples(object, use_names), Parameter %in% pars)
  if (object$coord_system != "proportions") {
    mtidy <- dplyr::group_by(mtidy, Parameter, coord, coord2, sample, covariate) 
  } else {
    mtidy <- dplyr::group_by(mtidy, Parameter, coord, sample, covariate)
  }
  if (!gather_prob){
    mtidy <- mtidy %>% 
      driver::summarise_posterior(val, ...) %>%
      dplyr::ungroup() %>%
      split(.$Parameter) %>% 
      purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
  } else if (gather_prob){
    mtidy <- mtidy %>% 
      dplyr::select(-iter) %>% 
      tidybayes::mean_qi(val, .prob=c(.5, .8, .95, .99)) %>% 
      dplyr::ungroup() %>% 
      split(.$Parameter) %>% 
      purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
  }
  return(mtidy)
}

#' Print dimensions and coordinate system information for mongrelfit object. 
#'
#' @param x an object of class mongrelfit
#' @param ... currently unused
#' @export
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' print(fit)
#' }
#' @seealso \code{\link{summary.mongrelfit}} summarizes posterior intervals 
print.mongrelfit <- function(x, ...){
  cat("Mongrelfit Object: \n" )
  cat(paste("  Number of Samples:\t\t", x$N, "\n"))
  cat(paste("  Number of Categories:\t\t", x$D, "\n"))
  cat(paste("  Number of Covariates:\t\t", x$Q, "\n"))
  cat(paste("  Number of Posterior Samples:\t", x$iter, "\n"))
  
  pars <- c("Eta", "Lambda", "Sigma")
  pars <- pars[pars %in% names(x)]
  pars <- paste(pars, collapse = "  ")
  cat(paste("  Contains Samples of Parameters:", pars, "\n"))
  
  if (x$coord_system=="alr"){
    cs <- x$alr_base
    nm <- x$names_categories
    if (!is.null(nm)) cs <- paste0(cs, " [", nm[x$alr_base], "]")
    cs <- paste("alr, reference category: ", cs)
  } else {
    cs <- x$coord_system
  }
  cat(paste("  Coordinate System:\t\t", cs))
}


#' Return regression coefficients of mongrelfit object
#' 
#' Returned as array of dimension (D-1) x Q x iter.
#' 
#' @param object an object of class mongrelfit
#' @param use_names if column and row names were passed for Y and X in 
#' call to \code{\link{mongrel}}, should these names be applied to output 
#' array. 
#' @return Array of dimension (D-1) x Q x iter
#' 
#' @export
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' coef(fit)
#' }
coef.mongrelfit <- function(object, use_names=TRUE){
  if (is.null(object$Lambda)) stop("mongrelfit object does not contain samples of Lambda")
  x <- object$Lambda
  if (use_names) return(name_array(x, object, list("cat", "cov", NULL)))
  return(x)
}

#' Convert object of class mongrelfit to a list
#' 
#' @param x an object of class mongrelfit
#' @param ... currently unused
#' 
#' @export
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' as.list(fit)
#' }
as.list.mongrelfit <- function(x,...){
  attr(x, "class") <- "list"
  return(x)
}

#' Predict response from new data
#' 
#' 
#' @param object An object of class mongrelfit
#' @param newdata An optional matrix for which to evaluate predictions. If NULL
#'   (default), the original data of the model is used. 
#' @param response Options = "LambdaX":Mean of regression, "Eta", "Y": counts
#' @param size the number of counts per sample if response="Y" (as vector or matrix), 
#'   default if newdata=NULL and response="Y" is to use colsums of m$Y. Otherwise
#'   uses median colsums of m$Y as default. If passed as a matrix should have dimensions
#'   ncol(newdata) x m$iter. 
#' @param use_names if TRUE apply names to output 
#' @param summary if TRUE, posterior summary of predictions are returned rather
#'   than samples
#' @param ... other arguments passed to summarise_posterior
#' 
#' @details currently only implmented for mongrelfit objects in coord_system "default"
#' "alr", or "ilr". 
#' 
#' @return (if summary==FALSE) array D x N x iter; (if summary==TRUE) 
#' tibble with calculated posterior summaries 
#' 
#' @export
#' @importFrom stats median predict runif
#' @examples 
#' sim <- mongrel_sim()
#' fit <- mongrel(sim$Y, sim$X)
#' predict(fit)
predict.mongrelfit <- function(object, newdata=NULL, response="LambdaX", size=NULL, 
                               use_names=TRUE, summary=FALSE, ...){
  
  if (!(object$coord_system %in% c("alr", "ilr"))){
    stop("currently only accepts mongrelfit objects in coord_system alr, or ilr")
  }

  if (is.null(newdata)) {
    newdata <- object$X
    if (response=="Y") size <-colSums(object$Y)
  } else {
    if ((response=="Y")&&(is.null(size))) size <- median(colSums(object$Y))
  }
  if ((response=="Y") && is.vector(size)){
    size <- replicate(object$iter, size)
  }
  
  # # Try to match rownames of newdata to avoid possible errors...
  # if (!is.null(rownames(newdata))) newdata <- newdata[object$names_covariates,]
  # would have error if names_covariates is NULL

  nnew <- ncol(newdata)
  
  # Draw LambdaX
  if (is.null(object$Lambda)) stop("mongrelfit object does not contain samples of Lambda")
  LambdaX <- array(0, dim = c(object$D-1, nnew, object$iter))
  for (i in 1:object$iter){
    LambdaX[,,i] <- object$Lambda[,,i] %*% newdata
  }
  if (use_names) LambdaX <- name_array(LambdaX, object,
                                              list("cat", colnames(newdata), 
                                                   NULL))
  if ((response == "LambdaX") && summary) {
    LambdaX <- gather_array(LambdaX, val, coord, sample, iter) %>% 
      group_by(coord, sample) %>% 
      summarise_posterior(val, ...) %>% 
      ungroup() %>% 
      name_tidy(object, list("coord" = "cat", "sample"=colnames(newdata)))
    return(LambdaX)
  }
  if (response == "LambdaX") return(LambdaX)
  
  # Draw Eta
  Eta <- array(0, dim=dim(LambdaX))
  zEta <- array(rnorm((object$D-1)*nnew*object$iter), dim = dim(Eta))
  for (i in 1:object$iter){
    Eta[,,i] <- LambdaX[,,i] + t(chol(object$Sigma[,,i]))%*%zEta[,,i]
  }
  if (use_names) Eta <- name_array(Eta, object, list("cat", colnames(newdata), 
                                                 NULL))
  if ((response=="Eta") && summary) {
    Eta <- gather_array(Eta, val, coord, sample, iter) %>% 
      group_by(coord, sample) %>% 
      summarise_posterior(val, ...) %>% 
      ungroup() %>% 
      name_tidy(object, list("coord" = "cat", "sample"=colnames(newdata)))
  }
  if (response=="Eta") return(Eta)
  
  # Draw Y
  if (is.null(object$Eta)) stop("mongrelfit object does not contain samples of Eta")
  
  com <- names(object)[!(names(object) %in% c("Lambda", "Sigma"))] # to save computation
  Pi <- mongrel_to_proportions(object[com])$Eta
  Ypred <- array(0, dim=c(object$D, nnew, object$iter))
  for (i in 1:object$iter){
    for (j in 1:nnew){
      Ypred[,j,i] <- rmultinom(1, size=size[j,i], prob=Pi[,j,i])
    }
  }
  if (use_names) name_array(Ypred, object, 
                                   list(object$names_categories, colnames(newdata), 
                                        NULL))
  if ((response == "Y") && summary) {
    Ypred <- gather_array(Ypred, val, coord, sample, iter) %>% 
      group_by(coord, sample) %>% 
      summarise_posterior(val, ...) %>% 
      ungroup() %>% 
      name_tidy(object, list("coord" = object$names_categories, 
                               "sample"= colnames(newdata)))
  }
  if (response=="Y") return(Ypred)
  stop("response parameter not recognized")
}


# access_dims -------------------------------------------------------------


#' @rdname access_dims
#' @export
ncategories.mongrelfit <- function(m){ m$D }

#' @rdname access_dims
#' @export
nsamples.mongrelfit <- function(m){ m$N }

#' @rdname access_dims
#' @export
ncovariates.mongrelfit <- function(m){ m$Q }



# name_dims ---------------------------------------------------------------

#' @rdname name_dims
#' @export
names_covariates.mongrelfit <- function(m){
  return(m$names_covariates)
}

#' @rdname name_dims
#' @export
names_samples.mongrelfit <- function(m){
  return(m$names_samples)
  
}


#' @rdname name_dims
#' @export
names_categories.mongrelfit <- function(m){
  return(m$names_categories)
  
}


#' @rdname name_dims
#' @export
`names_covariates<-.mongrelfit` <- function(m, value){
  if (!is.null(value)) stopifnot(m$Q == length(value))
  m$names_covariates <- value
  m <- name(m)
  return(m)
}

#' @rdname name_dims
#' @export
`names_samples<-.mongrelfit` <- function(m, value){
  if (!is.null(value)) stopifnot(m$N == length(value))
  m$names_samples <- value
  m <- name(m)
  return(m)
}


#' @rdname name_dims
#' @export
`names_categories<-.mongrelfit` <- function(m, value){
  if (!is.null(value)) stopifnot(m$D == length(value))
  m$names_categories <- value
  m <- name(m)
  return(m)
}


# sample_prior ------------------------------------------------------------

#' Sample from the prior distribution of mongrelfit object
#' 
#' Note this can be used to sample from prior and then predict can
#' be called to get counts or LambdaX (\code{\link{predict.mongrelfit}})
#' 
#' @param m object of class mongrelfit
#' @param n_sample number of samples to produce
#' @param pars parameters to sample
#' @param use_names should names be used if available
#' @param ... currently ignored
#' @export
#' @importFrom stats rWishart
#' 
#' @details Could be greatly speed up in the future if needed by sampling
#' directly from cholesky form of inverse wishart (currently implemented as 
#' header in this library - see MatDist.h).  
#' @examples 
#' # Sample prior of already fitted  mongrelfit object
#' sim <- mongrel_sim()
#' attach(sim)
#' fit <- mongrel(Y, X)
#' sample_prior(fit)
#' 
#' # Sample prior as part of model fitting
#' m <- mongrelfit(N=N, D=D, Q=Q, iter=2000, upsilon=upsilon, 
#'                 Xi=Xi, Gamma=Gamma, Theta=Theta, X=X, 
#'                 coord_system="alr", alr_base=D)
#' m <- sample_prior(mongrelfit)
#' plot(m) # plot prior distribution (defaults to parameter Lambda) 
sample_prior.mongrelfit <- function(m, n_sample=2000, 
                                    pars=c("Eta", "Lambda", "Sigma"), 
                                    use_names=TRUE, ...){
  req(m, c("upsilon", "Theta", "Gamma", "Xi"))
  
  # Convert to default ALR for computation
  l <- store_coord(m)
  m <- mongrel_to_alr(m, m$D)
  
  # Sample Priors - Sigma
  LSigmaInv <- rWishart(n_sample, m$upsilon, solve(m$Xi))
  for (i in 1:n_sample) LSigmaInv[,,i] <- t(chol(LSigmaInv[,,i]))
  
  # Sample Priors - Lambda
  if (any(c("Eta", "Lambda") %in% pars)){
    Lambda <- array(rnorm((m$D-1)*m$Q*n_sample), dim=c(m$D-1, m$Q, n_sample)) 
    UGamma <- chol(m$Gamma)
    for (i in 1:n_sample){
      Lambda[,,i] <- m$Theta + forwardsolve(LSigmaInv[,,i], Lambda[,,i]) %*% UGamma 
    }  
  }
  
  # Sample Priors - Eta 
  if ("Eta" %in% pars){
    req(m, "X")
    Eta <- array(rnorm((m$D-1)*m$N*n_sample), dim=c(m$D-1, m$N, n_sample))
    for (i in 1:n_sample) {
      Eta[,,i] <- Lambda[,,i] %*% m$X + forwardsolve(LSigmaInv[,,i], Eta[,,i])
    }
  }
  
  # Solve for Sigma if requested
  if ("Sigma" %in% pars){
    Sigma <- LSigmaInv # to make code more readable at memory expense
    for (i in 1:n_sample) Sigma[,,i] <- chol2inv(t(Sigma[,,i]))
  }
  
  # Convert to object of class mongrelfit
  out <- mongrelfit(m$D, m$N, m$Q, iter=as.integer(n_sample),
                    coord_system="alr", 
                    alr_base=m$D, 
                    Eta = mifelse("Eta" %in% pars, Eta, NULL), 
                    Sigma = mifelse("Sigma" %in% pars, Sigma, NULL), 
                    Lambda = mifelse("Lambda" %in% pars, Lambda, NULL), 
                    Xi = m$Xi, 
                    upsilon=m$upsilon, 
                    Theta = m$Theta, 
                    X = mifelse("Sigma" %in% pars, m$X, NULL), 
                    Gamma=m$Gamma, 
                    names_covariates=m$names_covariates, 
                    names_samples = m$names_samples, 
                    names_categories = m$names_categories)
  
  # Convert back to original afterwards
  out <- reapply_coord(out, l)
  if (use_names) out <- name(out)
  verify(out)
  return(out)
}
