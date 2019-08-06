#' Convert pibble samples of Eta Lambda and Sigma to tidy format
#' 
#' Combines them all into a single tibble, see example for formatting and 
#' column headers. Primarily designed to be used by 
#' \code{\link{summary.pibblefit}}. 
#' 
#' @param m an object of class pibblefit
#' @param use_names should dimension indices be replaced by
#'   dimension names if provided in data used to fit pibble model.  
#' @param as_factor if use_names should names be returned as factor?
#' 
#' @importFrom driver gather_array
#' @importFrom dplyr bind_rows group_by 
#' @export
#' @return tibble
#' @examples 
#' sim <- pibble_sim()
#' fit <- pibble(sim$Y, sim$X)
#' fit_tidy <- pibble_tidy_samples(fit, use_names=TRUE)
#' head(fit_tidy)
pibble_tidy_samples<- function(m, use_names=FALSE, as_factor=FALSE){
  l <- list()
  if (!is.null(m$Eta)) l$Eta <- driver::gather_array(m$Eta, .data$val, 
                                                     .data$coord, 
                                                     .data$sample, 
                                                     .data$iter)
  if (!is.null(m$Lambda)) l$Lambda <- driver::gather_array(m$Lambda, .data$val, 
                                                           .data$coord, 
                                                           .data$covariate, 
                                                           .data$iter)
  if (!is.null(m$Sigma)) l$Sigma <- driver::gather_array(m$Sigma, .data$val, 
                                                         .data$coord, 
                                                         .data$coord2, 
                                                         .data$iter) 

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

  l <- name_tidy(l, m, cl, as_factor)
  
  return(l)
}


#' Convert orthus samples of Eta Lambda and Sigma to tidy format
#' 
#' Combines them all into a single tibble, see example for formatting and 
#' column headers. Primarily designed to be used by 
#' \code{\link{summary.orthusfit}}. 
#' 
#' @param m an object of class orthusfit
#' @param use_names should dimension indices be replaced by
#'   dimension names if provided in data used to fit pibble model.  
#' @param as_factor if use_names should names be returned as factor?
#' 
#' @importFrom driver gather_array
#' @importFrom dplyr bind_rows group_by 
#' @export
#' @return tibble
#' @examples 
#' sim <- orthus_sim()
#' fit <- orthus(sim$Y, sim$Z, sim$X)
#' fit_tidy <- orthus_tidy_samples(fit, use_names=TRUE)
#' head(fit_tidy)
orthus_tidy_samples<- function(m, use_names=FALSE, as_factor=FALSE){
  l <- list()
  if (!is.null(m$Eta)) l$Eta <- driver::gather_array(m$Eta, .data$val, 
                                                     .data$coord, .data$sample, 
                                                     .data$iter)
  if (!is.null(m$Lambda)) l$Lambda <- driver::gather_array(m$Lambda, .data$val, 
                                                           .data$coord, .data$covariate, 
                                                           .data$iter)
  if (!is.null(m$Sigma)) l$Sigma <- driver::gather_array(m$Sigma, .data$val, 
                                                         .data$coord, .data$coord2, 
                                                         .data$iter) 
  
  l <- dplyr::bind_rows(l, .id="Parameter")
  
  if (!use_names) return(l)
  # Deal with names (first create conversion list)
  cl <- list()
  if (!is.null(m$Eta)) {
    cl[["sample"]] = "sam"
    cl[["coord"]] = "combo"
  }
  if (!is.null(m$Lambda)){
    cl[["covariate"]] = "cov"
    cl[["coord"]] = "combo"
  }
  if (!is.null(m$Sigma)){
    cl[["coord"]] = "combo"
    cl[["coord2"]] = "combo"
  }
  
  l <- name_tidy(l, m, cl, as_factor)
  
  return(l)
}



# Internal function to check if summary has already been precomputed. 
summary_check_precomputed <- function(m, pars){
  if (!is.null(m$summary)){
    if (all(!is.null(m$summary[pars]))) return(TRUE)
  }
  return(FALSE)
}

#' Summarise pibblefit object and print posterior quantiles
#' 
#' Default calculates median, mean, 50\% and 95\% credible interval
#' 
#' @param object an object of class pibblefit 
#' @param pars character vector (default: c("Eta", "Lambda", "Sigma"))
#' @param use_names should summary replace dimension indices with pibblefit 
#'   names if names Y and X were named in call to \code{\link{pibble}}
#' @param as_factor if use_names and as_factor then returns names as factors 
#'   (useful for maintaining orderings when plotting)
#' @param gather_prob if TRUE then prints quantiles in long format rather than 
#'  wide (useful for some plotting functions)
#' @param ... other expressions to pass to summarise (using name 'val' unquoted is 
#'   probably what you want)
#' @import dplyr
#' @importFrom driver summarise_posterior
#' @importFrom purrr map
#' @importFrom tidybayes mean_qi
#' @importFrom dplyr group_by select ungroup
#' @importFrom rlang syms
#' @export
#' @examples 
#' \dontrun{
#' fit <- pibble(Y, X)
#' summary(fit, pars="Eta", median = median(val))
#' 
#' # Some later functions make use of precomputation
#' fit$summary <- summary(fit)
#' }
summary.pibblefit <- function(object, pars=NULL, use_names=TRUE, as_factor=FALSE, 
                               gather_prob=FALSE, ...){
  if (is.null(pars)) {
    pars <- c()
    if (!is.null(object$Eta)) pars <- c(pars, "Eta")
    if (!is.null(object$Lambda)) pars <- c(pars, "Lambda")
    if (!is.null(object$Sigma)) pars <- c(pars, "Sigma")
    pars <- pars[pars %in% names(object)] # only for the ones that are present 
  }
  
  
  # if already calculated
  if (summary_check_precomputed(object, pars)) return(object$summary[pars])
  
  mtidy <- dplyr::filter(pibble_tidy_samples(object, use_names, as_factor), 
                         .data$Parameter %in% pars)
  # Suppress warnings about stupid implict NAs, this is on purpose. 
  suppressWarnings({
    
    vars <- c()
    if ("Eta" %in% pars) vars <- c(vars, "coord", "sample")
    if ("Lambda" %in% pars) vars <- c(vars, "coord", "covariate")
    if (("Sigma" %in% pars) & (object$coord_system != "proportions")) {
      vars <- c(vars, "coord", "coord2")
    }
    vars <- unique(vars)
    vars <- rlang::syms(vars)
    
    mtidy <- dplyr::group_by(mtidy, .data$Parameter, !!!vars)
    # if ((object$coord_system != "proportions")) {
    #   mtidy <- dplyr::group_by(mtidy, Parameter, coord, coord2, sample, covariate) 
    # } else {
    #   mtidy <- dplyr::group_by(mtidy, Parameter, coord, sample, covariate)
    # }
    if (!gather_prob){
      mtidy <- mtidy %>% 
        driver::summarise_posterior(.data$val, ...) %>%
        dplyr::ungroup() %>%
        split(.$Parameter) %>% 
        purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
    } else if (gather_prob){
      mtidy <- mtidy %>% 
        dplyr::select(-.data$iter) %>% 
        tidybayes::mean_qi(.data$val, .width=c(.5, .8, .95, .99)) %>% 
        dplyr::ungroup() %>% 
        split(.$Parameter) %>% 
        purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
    }
      
  })
  
  return(mtidy)
}

#' Summarise orthusfit object and print posterior quantiles
#' 
#' Default calculates median, mean, 50\% and 95\% credible interval
#' 
#' @param object an object of class orthusfit 
#' @param pars character vector (default: c("Eta", "Lambda", "Sigma"))
#' @param use_names should summary replace dimension indices with orthusfit 
#'   names if names Y and X were named in call to \code{\link{orthus}}
#' @param as_factor if use_names and as_factor then returns names as factors 
#'   (useful for maintaining orderings when plotting)
#' @param gather_prob if TRUE then prints quantiles in long format rather than 
#'  wide (useful for some plotting functions)
#' @param ... other expressions to pass to summarise (using name 'val' unquoted is 
#'   probably what you want)
#' @import dplyr
#' @importFrom driver summarise_posterior
#' @importFrom purrr map
#' @importFrom tidybayes mean_qi
#' @importFrom dplyr group_by select ungroup
#' @importFrom rlang syms
#' @export
#' @examples 
#' \dontrun{
#' fit <- orthus(Y, Z, X)
#' summary(fit, pars="Eta", median = median(val))
#' 
#' # Some later functions make use of precomputation
#' fit$summary <- summary(fit)
#' }
summary.orthusfit <- function(object, pars=NULL, use_names=TRUE, as_factor=FALSE, 
                              gather_prob=FALSE, ...){
  if (is.null(pars)) {
    pars <- c()
    if (!is.null(object$Eta)) pars <- c(pars, "Eta")
    if (!is.null(object$Lambda)) pars <- c(pars, "Lambda")
    if (!is.null(object$Sigma)) pars <- c(pars, "Sigma")
    pars <- pars[pars %in% names(object)] # only for the ones that are present 
  }
  
  
  # if already calculated
  if (summary_check_precomputed(object, pars)) return(object$summary[pars])
  
  mtidy <- dplyr::filter(orthus_tidy_samples(object, use_names, as_factor), 
                         .data$Parameter %in% pars)
  # Suppress warnings about stupid implict NAs, this is on purpose. 
  suppressWarnings({
    
    vars <- c()
    if ("Eta" %in% pars) vars <- c(vars, "coord", "sample")
    if ("Lambda" %in% pars) vars <- c(vars, "coord", "covariate")
    if (("Sigma" %in% pars) & (object$coord_system != "proportions")) {
      vars <- c(vars, "coord", "coord2")
    }
    vars <- unique(vars)
    vars <- rlang::syms(vars)
    
    mtidy <- dplyr::group_by(mtidy, .data$Parameter, !!!vars)
    # if ((object$coord_system != "proportions")) {
    #   mtidy <- dplyr::group_by(mtidy, Parameter, coord, coord2, sample, covariate) 
    # } else {
    #   mtidy <- dplyr::group_by(mtidy, Parameter, coord, sample, covariate)
    # }
    if (!gather_prob){
      mtidy <- mtidy %>% 
        driver::summarise_posterior(.data$val, ...) %>%
        dplyr::ungroup() %>%
        split(.$Parameter) %>% 
        purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
    } else if (gather_prob){
      mtidy <- mtidy %>% 
        dplyr::select(-.data$iter) %>% 
        tidybayes::mean_qi(.data$val, .width=c(.5, .8, .95, .99)) %>% 
        dplyr::ungroup() %>% 
        split(.$Parameter) %>% 
        purrr::map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
    }
    
  })
  
  return(mtidy)
}

#' Print dimensions and coordinate system information for pibblefit object. 
#'
#' @param x an object of class pibblefit
#' @param summary if true also calculates and prints summary
#' @param ... other arguments to pass to summary function
#' @export
#' @examples 
#' \dontrun{
#' fit <- pibble(Y, X)
#' print(fit)
#' }
#' @seealso \code{\link{summary.pibblefit}} summarizes posterior intervals 
print.pibblefit <- function(x, summary=FALSE, ...){
  if (is.null(x$Y)) {
    cat(" pibblefit Object (Priors Only): \n")
  } else {
    cat("pibblefit Object: \n" )  
  }
  
  cat(paste("  Number of Samples:\t\t", x$N, "\n"))
  cat(paste("  Number of Categories:\t\t", x$D, "\n"))
  cat(paste("  Number of Covariates:\t\t", x$Q, "\n"))
  cat(paste("  Number of Posterior Samples:\t", x$iter, "\n"))
  
  pars <- c("Eta", "Lambda", "Sigma")
  pars <- pars[pars %in% names(x)]
  pars <- paste(pars, collapse = "  ")
  cat(paste("  Contains Samples of Parameters:", pars, "\n", sep=""))
  
  if (x$coord_system=="alr"){
    cs <- x$alr_base
    nm <- x$names_categories
    if (!is.null(nm)) cs <- paste0(cs, " [", nm[x$alr_base], "]")
    cs <- paste("alr, reference category:", cs)
  } else {
    cs <- x$coord_system
  }
  cat(paste("  Coordinate System:\t\t", cs, "\n"))
  if (!is.null(x$logMarginalLikelihood)){
    cat(paste("  Log Marginal Likelihood:\t", 
              round(x$logMarginalLikelihood, 3), "\n"))
  }
  
  if (summary){
    cat("\n\n Summary: \n ")
    print(summary(x, ...))
  }
}


#' Print dimensions and coordinate system information for orthusfit object. 
#'
#' @param x an object of class orthusfit
#' @param summary if true also calculates and prints summary
#' @param ... other arguments to pass to summary function
#' @export
#' @examples 
#' \dontrun{
#' fit <- orthus(Y, Z, X)
#' print(fit)
#' }
#' @seealso \code{\link{summary.orthusfit}} summarizes posterior intervals 
print.orthusfit <- function(x, summary=FALSE, ...){
  if (is.null(x$Y)) {
    cat(" orthusfit Object (Priors Only): \n")
  } else {
    cat("orthusfit Object: \n" )  
  }
  
  cat(paste("  Number of Samples:\t\t", x$N, "\n"))
  cat(paste("  Number of Categories:\t\t", x$D, "\n"))
  cat(paste("  Number of Zdimensions:\t", x$P, "\n"))
  cat(paste("  Number of Covariates:\t\t", x$Q, "\n"))
  cat(paste("  Number of Posterior Samples:\t", x$iter, "\n"))
  
  pars <- c("Eta", "Lambda", "Sigma")
  pars <- pars[pars %in% names(x)]
  pars <- paste(pars, collapse = "  ")
  cat(paste("  Contains Samples of Parameters:", pars, "\n", sep=""))
  
  if (x$coord_system=="alr"){
    cs <- x$alr_base
    nm <- x$names_categories
    if (!is.null(nm)) cs <- paste0(cs, " [", nm[x$alr_base], "]")
    cs <- paste("alr, reference category:", cs)
  } else {
    cs <- x$coord_system
  }
  cat(paste("  Coordinate System:\t\t", cs, "\n"))
  if (!is.null(x$logMarginalLikelihood)){
    cat(paste("  Log Marginal Likelihood:\t", 
              round(x$logMarginalLikelihood, 3), "\n"))
  }
  
  if (summary){
    cat("\n\n Summary: \n ")
    print(summary(x, ...))
  }
}



#' Return regression coefficients of pibblefit object
#' 
#' Returned as array of dimension (D-1) x Q x iter (if in ALR or ILR) otherwise
#' DxQxiter (if in proportions or clr).
#' 
#' @param object an object of class pibblefit
#' @param ... other options passed to coef.pibblefit (see details)
#' @return Array of dimension (D-1) x Q x iter
#' @details Other arguments:
#' \itemize{
#' \item `use_names` if column and row names were passed for Y and X in 
#' call to \code{\link{pibble}}, should these names be applied to output 
#' array. 
#' }
#' 
#' @export
#' @examples 
#' \dontrun{
#' fit <- pibble(Y, X)
#' coef(fit)
#' }
coef.pibblefit <- function(object, ...){
  args <- list(...)
  use_names <- args_null("use_names", args, TRUE)
  
  if (is.null(object$Lambda)) stop("pibblefit object does not contain samples of Lambda")
  x <- object$Lambda
  if (use_names) return(name_array(x, object, list("cat", "cov", NULL)))
  return(x)
}


#' Return regression coefficients of orthus object
#' 
#' Returned as array of dimension (D-1+P) x Q x iter (if in ALR or ILR) 
#' otherwise (D+P) x Q x iter.
#' 
#' @param object an object of class orthusfit
#' @param ... other options passed to coef.orthusfit (see details)
#' @return Array of dimension (D-1) x Q x iter
#' @details Other arguments:
#' \itemize{
#' \item use_names if column and row names were passed for Y and X in 
#' call to \code{\link{pibble}}, should these names be applied to output 
#' array. 
#' }
#' 
#' @export
#' @examples 
#' \dontrun{
#' fit <- orthus(Y, Z, X)
#' coef(fit)
#' }
coef.orthusfit <- function(object, ...){
  args <- list(...)
  use_names <- args_null("use_names", args, TRUE)
  if (is.null(object$Lambda)) stop("orthusfit object does not contain samples of Lambda")
  if (use_names) object <- name.orthusfit(object)
  x <- object$Lambda
  return(x)
}


#' Convert object of class pibblefit to a list
#' 
#' @param x an object of class pibblefit
#' @param ... currently unused
#' 
#' @export
#' @examples 
#' \dontrun{
#' fit <- pibble(Y, X)
#' as.list(fit)
#' }
as.list.pibblefit <- function(x,...){
  attr(x, "class") <- "list"
  return(x)
}



#' Convert object of class orthusfit to a list
#' 
#' @param x an object of class orthusfit
#' @param ... currently unused
#' 
#' @export
#' @examples 
#' \dontrun{
#' fit <- orthus(Y, Z, X)
#' as.list(fit)
#' }
as.list.orthusfit <- function(x,...){
  attr(x, "class") <- "list"
  return(x)
}


#' Predict response from new data
#' 
#' 
#' @param object An object of class pibblefit
#' @param newdata An optional matrix for which to evaluate predictions. If NULL
#'   (default), the original data of the model is used. 
#' @param response Options = "LambdaX":Mean of regression, "Eta", "Y": counts
#' @param size the number of counts per sample if response="Y" (as vector or matrix), 
#'   default if newdata=NULL and response="Y" is to use colsums of m$Y. Otherwise
#'   uses median colsums of m$Y as default. If passed as a matrix should have dimensions
#'   ncol(newdata) x iter. 
#' @param use_names if TRUE apply names to output 
#' @param summary if TRUE, posterior summary of predictions are returned rather
#'   than samples
#' @param iter number of iterations to return if NULL uses object$iter
#' @param from_scratch should predictions of Y come from fitted Eta or from 
#'   predictions of Eta from posterior of Lambda? (default: false)
#' @param ... other arguments passed to summarise_posterior
#' 
#' @details currently only implemented for pibblefit objects in coord_system "default"
#' "alr", or "ilr". 
#' 
#' @return (if summary==FALSE) array D x N x iter; (if summary==TRUE) 
#' tibble with calculated posterior summaries 
#' 
#' @export
#' @importFrom stats median predict runif
#' @examples 
#' sim <- pibble_sim()
#' fit <- pibble(sim$Y, sim$X)
#' predict(fit)[,,1:2] # just show 2 samples
predict.pibblefit <- function(object, newdata=NULL, response="LambdaX", size=NULL, 
                               use_names=TRUE, summary=FALSE, iter=NULL, from_scratch=FALSE, ...){
  
  l <- store_coord(object)
  if (!(object$coord_system %in% c("alr", "ilr"))){
    object <- to_alr(object, ncategories(object))
    transformed <- TRUE
  } else {
    transformed <- FALSE
  }
  
  # If newdata is null - then predict based on existing data (X)
  # If size is null - then use colsums or median colsums of Y, 
  #    if that is null throw informative error
  if (is.null(newdata)){
    newdata <- object$X
    if (response=="Y"){
      if (is.null(size)){
        if (is.null(object$Y)){
          stop("Either Y or size must be specified to predict Counts")
        } else { # Y not null
          size <- colSums(object$Y)
        }
      }
    }
  } else { #newdata specified 
    if (response=="Y"){
      if (is.null(size)){
        if (is.null(object$Y)){
          stop("Either Y or size must be specified to predict Counts")
        } else {
          size <- median(colSums(object$Y))
        }
      }
    }
  }
  
  # if iter is null use object$iter
  if (is.null(iter)){ iter <- object$iter }
  
  # if size is a scalar, replicate it to a vector 
  if ((response=="Y") && (length(size)==1)) { size <- replicate(ncol(newdata), size) }
  # If size is a vector, replicate it to a matrix
  if ((response=="Y") && is.vector(size)){ size <- replicate(iter, size) }

  
  nnew <- ncol(newdata)
  
  # Draw LambdaX
  if (is.null(object$Lambda)) stop("pibblefit object does not contain samples of Lambda")
  LambdaX <- array(0, dim = c(object$D-1, nnew, iter))
  for (i in 1:iter){
    LambdaX[,,i] <- object$Lambda[,,i] %*% newdata
  }
  if (use_names) LambdaX <- name_array(LambdaX, object,
                                       list("cat", colnames(newdata), 
                                            NULL))
  if (response=="LambdaX"){
    if (transformed){
      LambdaX <- alrInv_array(LambdaX, object$D, 1)
      if (l$coord_system == "clr") LambdaX <- clr_array(LambdaX, 1)
    }
  }
  if ((response == "LambdaX") && summary) {
    LambdaX <- gather_array(LambdaX, .data$val, .data$coord, .data$sample, .data$iter) %>% 
      group_by(.data$coord, .data$sample) %>% 
      summarise_posterior(.data$val, ...) %>% 
      ungroup() %>% 
      name_tidy(object, list("coord" = "cat", "sample"=colnames(newdata)))
    return(LambdaX)
  }
  if (response == "LambdaX") return(LambdaX)
  
  # Draw Eta
  Eta <- array(0, dim=dim(LambdaX))
  zEta <- array(rnorm((object$D-1)*nnew*iter), dim = dim(Eta))
  for (i in 1:iter){
    Eta[,,i] <- LambdaX[,,i] + t(chol(object$Sigma[,,i]))%*%zEta[,,i]
  }
  if (use_names) Eta <- name_array(Eta, object, list("cat", colnames(newdata), 
                                                     NULL))
  if (response=="Eta"){
    if (transformed){
      Eta <- alrInv_array(Eta, object$D, 1)
      if (l$coord_system == "clr") Eta <- clr_array(Eta, 1)
    }
  }
  if ((response=="Eta") && summary) {
    Eta <- gather_array(Eta, .data$val, .data$coord, .data$sample, .data$iter) %>% 
      group_by(.data$coord, .data$sample) %>% 
      summarise_posterior(.data$val, ...) %>% 
      ungroup() %>% 
      name_tidy(object, list("coord" = "cat", "sample"=colnames(newdata)))
  }
  if (response=="Eta") return(Eta)
  
  # Draw Y
  if (from_scratch){
    Pi <- alrInv_array(Eta, d=nrow(Eta)+1, coords=1)
  } else {
    if (is.null(object$Eta)) stop("pibblefit object does not contain samples of Eta")
    com <- names(object)[!(names(object) %in% c("Lambda", "Sigma"))] # to save computation
    Pi <- to_proportions(as.pibblefit(object[com]))$Eta
  }
  Ypred <- array(0, dim=c(object$D, nnew, iter))
  for (i in 1:iter){
    for (j in 1:nnew){
      Ypred[,j,i] <- rmultinom(1, size=size[j,i], prob=Pi[,j,i])
    }
  }
  if (use_names) name_array(Ypred, object, 
                            list(object$names_categories, colnames(newdata), 
                                 NULL))
  if ((response == "Y") && summary) {
    Ypred <- gather_array(Ypred, .data$val, .data$coord, .data$sample, .data$iter) %>% 
      group_by(.data$coord, .data$sample) %>% 
      summarise_posterior(.data$val, ...) %>% 
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
ncategories.pibblefit <- function(m){ m$D }

#' @rdname access_dims
#' @export
nsamples.pibblefit <- function(m){ m$N }

#' @rdname access_dims
#' @export
ncovariates.pibblefit <- function(m){ m$Q }

#' @rdname access_dims
#' @export
niter.pibblefit <- function(m){ m$iter }


#' @rdname access_dims
#' @export
ncategories.orthusfit <- function(m){ m$D }

#' @rdname access_dims
#' @export
nsamples.orthusfit <- function(m){ m$N }

#' @rdname access_dims
#' @export
ncovariates.orthusfit <- function(m){ m$Q }

#' @rdname access_dims
#' @export
niter.orthusfit <- function(m){ m$iter }



# name_dims ---------------------------------------------------------------

#' @rdname name_dims
#' @export
names_covariates.pibblefit <- function(m){
  return(m$names_covariates)
}

#' @rdname name_dims
#' @export
names_samples.pibblefit <- function(m){
  return(m$names_samples)
  
}


#' @rdname name_dims
#' @export
names_categories.pibblefit <- function(m){
  return(m$names_categories)
  
}

#' @rdname name_dims
#' @export
names_coords.pibblefit <- function(m){
  return(assign_cat_names(m))
}



#' @rdname name_dims
#' @export
`names_covariates<-.pibblefit` <- function(m, value){
  if (!is.null(value)) stopifnot(m$Q == length(value))
  m$names_covariates <- value
  m <- name(m)
  return(m)
}

#' @rdname name_dims
#' @export
`names_samples<-.pibblefit` <- function(m, value){
  if (!is.null(value)) stopifnot(m$N == length(value))
  m$names_samples <- value
  m <- name(m)
  return(m)
}


#' @rdname name_dims
#' @export
`names_categories<-.pibblefit` <- function(m, value){
  if (!is.null(value)) stopifnot(m$D == length(value))
  m$names_categories <- value
  m <- name(m)
  return(m)
}


# sample_prior ------------------------------------------------------------

#' Sample from the prior distribution of pibblefit object
#' 
#' Note this can be used to sample from prior and then predict can
#' be called to get counts or LambdaX (\code{\link{predict.pibblefit}})
#' 
#' @param m object of class pibblefit
#' @param n_samples number of samples to produce
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
#' # Sample prior of already fitted  pibblefit object
#' sim <- pibble_sim()
#' attach(sim)
#' fit <- pibble(Y, X)
#' sample_prior(fit)
#' 
#' # Sample prior as part of model fitting
#' m <- pibblefit(N=as.integer(sim$N), D=as.integer(sim$D), Q=as.integer(sim$Q), 
#'                 iter=2000L, upsilon=upsilon, 
#'                 Xi=Xi, Gamma=Gamma, Theta=Theta, X=X, 
#'                 coord_system="alr", alr_base=D)
#' m <- sample_prior(m)
#' plot(m) # plot prior distribution (defaults to parameter Lambda) 
sample_prior.pibblefit <- function(m, n_samples=2000L, 
                                    pars=c("Eta", "Lambda", "Sigma"), 
                                    use_names=TRUE, ...){
  req(m, c("upsilon", "Theta", "Gamma", "Xi"))
  
  # Convert to default ALR for computation
  l <- store_coord(m)
  m <- to_alr(m, m$D)
  
  # Sample Priors - Sigma
  USigmaInv <- rWishart(n_samples, m$upsilon, solve(m$Xi))
  for (i in 1:n_samples) USigmaInv[,,i] <- chol(USigmaInv[,,i])
  
  # Sample Priors - Lambda
  if (any(c("Eta", "Lambda") %in% pars)){
    Lambda <- array(rnorm((m$D-1)*m$Q*n_samples), dim=c(m$D-1, m$Q, n_samples)) 
    UGamma <- chol(m$Gamma)
    for (i in 1:n_samples){
      Lambda[,,i] <- m$Theta + backsolve(USigmaInv[,,i], Lambda[,,i]) %*% UGamma 
    }  
  }
  
  # Sample Priors - Eta 
  if ("Eta" %in% pars){
    req(m, "X")
    Eta <- array(rnorm((m$D-1)*m$N*n_samples), dim=c(m$D-1, m$N, n_samples))
    for (i in 1:n_samples) {
      Eta[,,i] <- Lambda[,,i] %*% m$X + backsolve(USigmaInv[,,i], Eta[,,i])
    }
  }
  
  # Solve for Sigma if requested
  if ("Sigma" %in% pars){
    Sigma <- USigmaInv # to make code more readable at memory expense
    for (i in 1:n_samples) {
      Sigma[,,i] <- backsolve(Sigma[,,i], diag(m$D-1))
      Sigma[,,i] <- tcrossprod(Sigma[,,i])
    }
  }
  
  # Convert to object of class pibblefit
  out <- pibblefit(m$D, m$N, m$Q, iter=as.integer(n_samples),
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



# ppc_summary -------------------------------------------------------------

#' @rdname ppc_summary
#' @param from_scratch should predictions of Y come from fitted Eta or from 
#'   predictions of Eta from posterior of Lambda? (default: false)
#' @export
ppc_summary.pibblefit <- function(m, from_scratch=FALSE, ...){
  if (!is.null(m$Y)) {
    o <- order(m$Y, decreasing=TRUE)
  } else {
    stop("ppc_summary is only for posterior samples, current object has Y==NULL")
  }
  
  if (m$iter ==1){
    warning("ppc_summary is intended to be used with more than 1 summary, ", 
            "results will be missleading")
  }
  
  pp <- predict(m, response="Y", from_scratch=from_scratch)
  pp <- matrix(pp, m$D*m$N, m$iter) 
  pp <- pp[o,]
  
  tr <- data.frame(dim_1 = 1:(m$N*m$D), 
                   dim_2 = NA, 
                   val = c(m$Y)[o])

  pp <- apply(pp, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  rownames(pp) <- c("p2.5", "p97.5")
  pp <- as.data.frame(t(pp))
  pp$dim_1 <- 1:nrow(pp)
  
  inBounds <- (pp$p2.5 <= tr$val) & (pp$p97.5 >= tr$val)
  inBounds <- sum(inBounds)/length(inBounds)
  cat("Proportions of Observations within 95% Credible Interval: ")
  cat(inBounds)
  cat("\n")
  invisible(c("percent.in.p95"=inBounds))
}

