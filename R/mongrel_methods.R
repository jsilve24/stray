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

  l <- apply_names_tidy(l, m, cl)
  
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
#' Default calculates median, mean, 50% and 95% credible interval
#' 
#' @param m an object of class mongrelfit 
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
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' summary(fit, pars="Eta", median = median(val))
#' 
#' # Some later functions make use of precomputation
#' fit$summary <- summary(fit)
#' }
summary.mongrelfit <- function(m, pars=NULL, use_names=TRUE, gather_prob=FALSE,
                               ...){
  if (is.null(pars)) {
    pars <- c("Eta", "Lambda", "Sigma")
    pars <- pars[pars %in% names(m)] # only for the ones that are present 
  }
  
  
  # if already calculated
  if (summary_check_precomputed(m, pars)) return(m$summary[pars])
  
  mtidy <- dplyr::filter(mongrel_tidy_samples(m, use_names), Parameter %in% pars)
  if (m$coord_system != "proportions") {
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
      tidybayes::mean_qi(.prob=c(.5, .8, .95, .99)) %>% 
      dplyr::ungroup() %>% 
      split(.$Parameter) %>% 
      purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
  }
  return(mtidy)
}

#' Print dimensions and coordinate system information for mongrelfit object. 
#'
#' @param m an object of class mongrelfit
#' @export
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' print(fit)
#' }
#' @seealso \code{\link{summary.mongrelfit}} summarizes posterior intervals 
print.mongrelfit <- function(m){
  cat("Mongrelfit Object: \n" )
  cat(paste("  Number of Samples:\t\t", m$N, "\n"))
  cat(paste("  Number of Categories:\t\t", m$D, "\n"))
  cat(paste("  Number of Covariates:\t\t", m$Q, "\n"))
  cat(paste("  Number of Posterior Samples:\t", m$iter, "\n"))
  
  pars <- c("Eta", "Lambda", "Sigma")
  pars <- pars[pars %in% names(m)]
  pars <- paste(pars, collapse = "  ")
  cat(paste("  Contains Samples of Parameters:", pars, "\n"))
  
  if (m$coord_system=="alr"){
    cs <- m$alr_base
    nm <- m$names_categories
    if (!is.null(nm)) cs <- paste0(cs, " [", nm[m$alr_base], "]")
    cs <- paste("alr, reference category: ", cs)
  } else {
    cs <- m$coord_system
  }
  cat(paste("  Coordinate System:\t\t", cs))
}


#' Return regression coefficients of mongrelfit object
#' 
#' Returned as array of dimension (D-1) x Q x iter.
#' 
#' @param m an object of class mongrelfit
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
coef.mongrelfit <- function(m, use_names=TRUE){
  if (is.null(m$Lambda)) stop("mongrelfit object does not contain samples of Lambda")
  x <- m$Lambda
  if (use_names) return(apply_names_array(x, m, list("cat", "cov", NULL)))
  return(x)
}

#' Convert object of class mongrelfit to a list
#' 
#' @param m an object of class mongrelfit
#' 
#' @export
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' as.list(fit)
#' }
as.list.mongrelfit <- function(m){
  attr(m, "class") <- "list"
  return(m)
}

#' Predict response from new data
#' 
#' 
#' @param m An object of class mongrelfit
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
#' @export
#' @importFrom stats median predict runif
#' @examples 
#' sim <- mongrel_sim()
#' fit <- mongrel(sim$Y, sim$X)
#' predict(fit)
predict.mongrelfit <- function(m, newdata=NULL, response="LambdaX", size=NULL, 
                               use_names=TRUE, summary=FALSE, ...){
  
  if (!(m$coord_system %in% c("alr", "ilr"))){
    stop("currently only accepts mongrelfit objects in coord_system alr, or ilr")
  }

  if (is.null(newdata)) {
    newdata <- m$X
    if (response=="Y") size <-colSums(m$Y)
  } else {
    if ((response=="Y")&&(is.null(size))) size <- median(colSums(m$Y))
  }
  if ((response=="Y") && is.vector(size)){
    size <- replicate(m$iter, size)
  }
  
  # Try to match rownames of newdata to avoid possible errors...
  if (!is.null(rownames(newdata))) newdata <- newdata[m$names_covariates,]
  
  nnew <- ncol(newdata)
  
  # Draw LambdaX
  if (is.null(m$Lambda)) stop("mongrelfit object does not contain samples of Lambda")
  LambdaX <- array(0, dim = c(m$D-1, nnew, m$iter))
  for (i in 1:m$iter){
    LambdaX[,,i] <- m$Lambda[,,i] %*% newdata
  }
  if (use_names) LambdaX <- apply_names_array(LambdaX, m,
                                              list("cat", colnames(newdata), 
                                                   NULL))
  if ((response == "LambdaX") && summary) {
    LambdaX <- gather_array(LambdaX, val, coord, sample, iter) %>% 
      group_by(coord, sample) %>% 
      summarise_posterior(val, ...) %>% 
      ungroup() %>% 
      apply_names_tidy(m, list("coord" = "cat", "sample"=colnames(newdata)))
    return(LambdaX)
  }
  if (response == "LambdaX") return(LambdaX)
  
  # Draw Eta
  Eta <- array(0, dim=dim(LambdaX))
  zEta <- array(rnorm((m$D-1)*nnew*m$iter), dim = dim(Eta))
  for (i in 1:m$iter){
    Eta[,,i] <- LambdaX[,,i] + t(chol(m$Sigma[,,i]))%*%zEta[,,i]
  }
  if (use_names) Eta <- apply_names_array(Eta, m, list("cat", colnames(newdata), 
                                                 NULL))
  if ((response=="Eta") && summary) {
    Eta <- gather_array(Eta, val, coord, sample, iter) %>% 
      group_by(coord, sample) %>% 
      summarise_posterior(val, ...) %>% 
      ungroup() %>% 
      apply_names_tidy(m, list("coord" = "cat", "sample"=colnames(newdata)))
  }
  if (response=="Eta") return(Eta)
  
  # Draw Y
  if (is.null(m$Eta)) stop("mongrelfit object does not contain samples of Eta")
  
  com <- names(m)[!(names(m) %in% c("Lambda", "Sigma"))] # to save computation
  Pi <- mongrel_to_proportions(m[com])$Eta
  Ypred <- array(0, dim=c(m$D, nnew, m$iter))
  for (i in 1:m$iter){
    for (j in 1:nnew){
      Ypred[,j,i] <- rmultinom(1, size=size[j,i], prob=Pi[,j,i])
    }
  }
  if (use_names) apply_names_array(Ypred, m, 
                                   list(m$names_categories, colnames(newdata), 
                                        NULL))
  if ((response == "Y") && summary) {
    Ypred <- gather_array(Ypred, val, coord, sample, iter) %>% 
      group_by(coord, sample) %>% 
      summarise_posterior(val, ...) %>% 
      ungroup() %>% 
      apply_names_tidy(m, list("coord" = m$names_categories, 
                               "sample"= colnames(newdata)))
  }
  if (response=="Y") return(Ypred)
  stop("response parameter not recognized")
}


#' Simple helper functions to access mongrel fit dimensions
#'
#' @param m An object of class mongrelfit 
#' @details An alternative approach to accessing these dimensions is to 
#'   access them directly from the mongrelfit object using list indexing. 
#' * \code{ncategories} is equivalent to \code{m$D}
#' * \code{nsamples} is equivalent to \code{m$N}
#' * \code{ncovariates} is equivalent to \code{m$Q}
#' @return integer 
#' @name access_dims
#' @examples 
#' \dontrun{
#' m <- mongrel(Y, X)
#' ncategories(m)
#' nsamples(m)
#' ncovariates(m)
#' }
NULL

#' @rdname access_dims
#' @export
ncategories <- function(m){ m$D }

#' @rdname access_dims
#' @export
nsamples <- function(m){ m$N }

#' @rdname access_dims
#' @export
ncovariates <- function(m){ m$Q }


