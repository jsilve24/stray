#' Simple verification of passed bassetfit object
#' @param m an object of class bassetfit
#' @param ... not used
#' @return throws error if any verification tests fail
#' @export 
verify.bassetfit <- function(m, ...){
  # check basic dimensions that must always be present
  stopifnot(is.integer(m$N), is.integer(m$Q),
            is.integer(m$D))
  stopifnot(is.character(m$coord_system))
  
  if (m$coord_system == "ilr") stopifnot(!is.null(m$ilr_base))
  if (m$coord_system == "alr") stopifnot(!is.null(m$alr_base))
  if (m$coord_system != "proportions") stopifnot(is.null(m$Sigma_default))
  if (m$coord_system != "proportions") stopifnot(is.null(m$Xi_default))
  
  N <- m$N; D <- m$D; Q <- m$Q; iter <- m$iter
  Dm1 <- ifelse (m$coord_system %in% c("ilr", "alr"), D-1, D)
  
  # throw error if iter is null but Eta, Sigma, and Lambda are not. 
  if (is.null(m$iter)) stopifnot(is.null(m$Eta), is.null(m$Lambda), 
                                 is.null(m$Sigma), is.null(m$Sigma_default))
  ifnotnull(m$iter, stopifnot(is.integer(m$iter)))
  ifnotnull(m$Eta,check_dims(m$Eta, c(Dm1, N, iter), "bassetfit param Eta"))
  ifnotnull(m$Lambda,check_dims(m$Lambda, c(Dm1, N, iter), "bassetfit param Lambda"))
  ifnotnull(m$Sigma,check_dims(m$Sigma, c(Dm1, Dm1, iter), "bassetfit param Sigma"))
  ifnotnull(m$Sigma_default,check_dims(m$Sigma_default, c(D-1, D-1, iter), "bassetfit param Sigma_default"))
  ifnotnull(m$Y,check_dims(m$Y, c(D, N), "bassetfit param Y"))
  ifnotnull(m$X,check_dims(m$X, c(Q, N), "bassetfit param X"))
  ifnotnull(m$upsilon,check_dims(m$upsilon, c(1), "bassetfit param upsilon"))
  ifnotnull(m$Xi,check_dims(m$Xi, c(Dm1,Dm1), "bassetfit param Xi"))
  ifnotnull(m$Xi_default,check_dims(m$Xi_default, c(D-1, D-1), "bassetfit param Xi_default"))
  ifnotnull(m$init,check_dims(m$init, c(Dm1, N), "bassetfit param init"))
  ifnotnull(m$names_categories,check_dims(m$names_categories, c(D), "bassetfit param names_categories"))
  ifnotnull(m$names_samples,check_dims(m$names_samples, c(N), "bassetfit param names_samples"))
  ifnotnull(m$names_covariates,check_dims(m$names_covariates, c(Q), "bassetfit param names_covariates"))
}


#' Predict using basset
#' 
#' @param object An object of class pibblefit
#' @param newdata An optional matrix for which to evaluate prediction. 
#' @param response Options = "Lambda":Mean of regression, "Eta", "Y": counts
#' @param size the number of counts per sample if response="Y" (as vector or matrix), 
#'   default if newdata=NULL and response="Y" is to use colsums of m$Y. Otherwise
#'   uses median colsums of object$Y as default. If passed as a matrix should have dimensions
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
predict.bassetfit <- function(object, newdata, response="Lambda", size=NULL, 
                              use_names=TRUE, summary=FALSE, iter=NULL,
                              from_scratch=FALSE, ...){
  req(object, c("Lambda", "Sigma"))
  newdata <- vec_to_mat(newdata)
  
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
  if (response=="Y"){
    if (is.null(size)){
      if (is.null(object$Y)){
        stop("Either Y or size must be specified to predict Counts")
      } else {
        size <- median(colSums(object$Y))
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
  
  # Set up Function Evaluation
  obs <- c(rep(TRUE, ncol(object$X)), rep(FALSE, nnew)) 
  Theta <- object$Theta(cbind(object$X, newdata))
  Gamma <- object$Gamma(cbind(object$X, newdata))
  
  # Predict Lambda
  Gamma_oo <- Gamma[obs, obs, drop=F]
  Gamma_ou <- Gamma[obs, !obs, drop=F]
  Gamma_uu <- Gamma[!obs, !obs, drop=F]
  Gamma_ooIou <- solve(Gamma_oo, Gamma_ou)
  Gamma_schur <- Gamma_uu - t(Gamma_ou) %*% Gamma_ooIou 
  U_Gamma_schur <- chol(Gamma_schur)
  Theta_o <- Theta[,obs, drop=F]
  Theta_u <- Theta[,!obs, drop=F]
  
  Lambda_u <- array(0, dim=c(object$D-1, nnew, object$iter))
  
  # function for prediction - sample one Lambda_u
  lu <- function(Lambda_o, Sigma){
    Z <- matrix(rnorm(nrow(Gamma_uu)*(object$D-1)), object$D-1, nrow(Gamma_uu))
    Theta_u + (Lambda_o-Theta_o)%*%Gamma_ooIou + t(chol(Sigma))%*%Z%*%U_Gamma_schur
  }
  
  # Fill in and predict
  for (i in 1:object$iter){
    Lambda_u[,,i] <- lu(object$Lambda[,,i], object$Sigma[,,i])
  }
  
  if (use_names) Lambda_u <- name_array(Lambda_u, object,
                                       list("cat", colnames(newdata),NULL))
  
  if (response=="Lambda"){
    if (transformed){
      Lambda_u <- alrInv_array(Lambda_u, object$D, 1)
      if (l$coord_system == "clr") Lambda_u <- clr_array(Lambda_u, 1)
    }
  }
  if ((response == "Lambda") && summary) {
    Lambda_u <- gather_array(Lambda_u, .data$val, .data$coord, .data$sample, .data$iter) %>% 
      group_by(.data$coord, .data$sample) %>% 
      summarise_posterior(.data$val, ...) %>% 
      ungroup() %>% 
      name_tidy(object, list("coord" = "cat", "sample"=colnames(newdata)))
    return(Lambda_u)
  }
  if (response == "Lambda") return(Lambda_u)
  
  # Draw Eta
  Eta <- array(0, dim=dim(Lambda_u))
  zEta <- array(rnorm((object$D-1)*nnew*iter), dim = dim(Eta))
  for (i in 1:iter){
    Eta[,,i] <- Lambda_u[,,i] + t(chol(object$Sigma[,,i]))%*%zEta[,,i]
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
  if (!from_scratch){
    Pi <- alrInv_array(Eta, d=nrow(Eta)+1, coords=1)
  } else {
    if (is.null(object$Eta)) stop("pibblefit object does not contain samples of Eta")
    
    com <- names(object)[!(names(object) %in% c("Lambda", "Sigma"))] # to save computation
    Pi <- to_proportions(object[com])$Eta
    Pi <- alrInv_array(Eta, d=nrow(Eta)+1, coords=1)
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