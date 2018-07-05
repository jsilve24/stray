# convert mongrel samples of Eta Lambda and Sigma to tidy format
#' @importFrom driver gather_array
#' @importFrom dplyr bind_rows
mongrel_tidy_samples<- function(m, use_names=FALSE){
  l <- list()
  l$Eta <- gather_array(m$Eta, val, coord, sample, iter)
  l$Lambda <- gather_array(m$Lambda, val, coord, covariate, iter)
  if (m$coord_system != "proportions"){
    l$Sigma <- gather_array(m$Sigma, val, coord, coord2, iter) 
  } 
  l <- dplyr::bind_rows(l, .id="Parameter")
  
  # a bunch of crazziness to deal with possible names given 
  l <- apply_names_tidy(l, m, list("sample" = "sam", 
                                   "covariate" = "cov", 
                                   "coord" = "cat"))
  if (m$coord_system != "proportions"){
    l <- apply_names_tidy(l, m, list("coord2" = "cat"))
  }
  return(l)
}

summary_check_precomputed <- function(m, pars){
  if (!is.null(m$summary)){
    if (all(!is.null(m$summary[pars]))) return(TRUE)
  }
  return(FALSE)
}


#' @param ... other expressions to pass to summarise (using name 'val' unquoted is 
#'   probably what you want)
#' @import dplyr
#' @importFrom driver summarise_posterior
#' @importFrom purrr map
#' @importFrom tidybayes mean_qi
#' @examples 
#' ... what is fit ... # TODO
#' summary(fit, pars="Eta", median = median(val))
#' ... note allows precomputation as 
#' fit$summary <- summary(fit) to save time for future functions. 
summary.mongrelfit <- function(m, pars=NULL, use_names=TRUE, gather_prob=FALSE,
                               ...){
  if (is.null(pars)) pars <- c("Eta", "Lambda", "Sigma")
  
  # if already calculated
  if (summary_check_precomputed(m, pars)) return(m$summary[pars])
  
  mtidy <- dplyr::filter(mongrel_tidy_samples(m, use_names), Parameter %in% pars)
  if (m$coord_system != "proportions") {
    mtidy <- group_by(mtidy, Parameter, coord, coord2, sample, covariate) 
  } else {
    mtidy <- group_by(mtidy, Parameter, coord, sample, covariate)
  }
  if (!gather_prob){
    mtidy <- mtidy %>% 
      summarise_posterior(val, ...) %>%
      ungroup() %>%
      split(.$Parameter) %>% 
      map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
  } else if (gather_prob){
    mtidy <- mtidy %>% 
      select(-iter) %>% 
      tidybayes::mean_qi(.prob=c(.5, .8, .95, .99)) %>% 
      ungroup() %>% 
      split(.$Parameter) %>% 
      map(~dplyr::select_if(.x, ~!all(is.na(.x))))  
  }
  return(mtidy)
}


print.mongrelfit <- function(m){
  cat("Mongrelfit Object: \n" )
  cat(paste("  Number of Samples:\t\t", m$N, "\n"))
  cat(paste("  Number of Categories:\t\t", m$D, "\n"))
  cat(paste("  Number of Covariates:\t\t", m$Q, "\n"))
  cat(paste("  Number of Posterior Samples:\t", m$iter, "\n"))
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


coef.mongrelfit <- function(m, use.names=TRUE){
  x <- m$Lambda
  if (use.names) return(apply_names_array(x, m, list("cat", "cov", NULL)))
  return(x)
}

coefficients.mongrelfit <- function(m, use_names=TRUE){
  coef.mongrelfit(m, use_names=TRUE)
}

as.list.mongrelfit <- function(m){
  attr(m, "class") <- "list"
  return(m)
}

#' Predict response from new data
#' 
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
#' 
predict.mongrelfit <- function(m, newdata=NULL, response="LambdaX", size=NULL, 
                               use_names=TRUE, summary=FALSE, ...){
  if (!(m$coord_system %in% c("alr", "ilr"))){
    stop("currently only accepts mongrelfit objects in coord_system default, alr, or ilr")
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
  Pi <- mongrel_to_proportions(m)$Eta
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


