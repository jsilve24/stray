#' Interface to fit orthus models 
#' 
#' This function is largely a more user friendly wrapper around 
#' \code{\link{optimPibbleCollapsed}} and 
#' \code{\link{uncollapsePibble}} for fitting orthus models. 
#' See details for model specification. 
#'  Notation: \code{N} is number of samples, \code{P} is the number of dimensions
#'  of observations in the second dataset,
#'  \code{D} is number of multinomial categories, \code{Q} is number
#'  of covariates, \code{iter} is the number of samples of \code{eta} (e.g.,
#'  the parameter \code{n_samples} in the function 
#'  \code{\link{optimPibbleCollapsed}})
#' @param Y D x N matrix of counts (if NULL uses priors only)
#' @param Z P x N matrix of counts (if NULL uses priors only - must be present/absent 
#'  if Y is present/absent)
#' @param X Q x N matrix of covariates (design matrix) (if NULL uses priors only, must
#' be present to sample Eta)
#' @param upsilon dof for inverse wishart prior (numeric must be > D) 
#'   (default: D+3)
#' @param Theta (D-1+P) x Q matrix of prior mean for regression parameters
#'   (default: matrix(0, D-1+P, Q))
#' @param Gamma QxQ prior covariance matrix 
#'   (default: diag(Q))
#' @param Xi (D-1+P)x(D-1+P) prior covariance matrix
#'   (default: ALR transform of diag(1)*(upsilon-D)/2 - this is 
#'   essentially iid on "base scale" using Aitchison terminology)
#' @param init (D-1) x Q initialization for Eta for optimization
#' @param pars character vector of posterior parameters to return
#' @param ... arguments passed to \code{\link{optimPibbleCollapsed}} and 
#'   \code{\link{uncollapsePibble}}
#' 
#' @details the full model is given by:
#'    \deqn{Y_j \sim Multinomial(Pi_j)} 
#'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
#'    \deqn{cbind(Eta, Z) \sim MN_{D-1+P x N}(Lambda*X, Sigma, I_N)}
#'    \deqn{Lambda \sim MN_{D-1+P x Q}(Theta, Sigma, Gamma)}
#'    \deqn{Sigma \sim InvWish(upsilon, Xi)}
#'  Where Gamma is a Q x Q covariance matrix, and Phi^{-1} is 
#'  ALRInv_D transform. 
#'  That is, the orthus model models the latent multinomial log-ratios (Eta) and
#'  the observations of the second dataset jointly as a linear model. This allows 
#'  Sigma to also describe the covariation between the two datasets.  
#'  
#'  Default behavior is to use MAP estimate for uncollaping the LTP 
#'  model if laplace approximation is not preformed. 
#' @return an object of class pibblefit
#' @md
#' @name orthus_fit
#' @examples 
#' sim <- orthus_sim()
#' fit <- orthus(sim$Y, sim$Z, sim$X)
#' @seealso \code{\link{stray_transforms}} provide convenience methods for 
#'  transforming the representation of pibblefit objects (e.g., conversion to 
#'  proportions, alr, clr, or ilr coordinates.)
#'  
#' \code{\link{access_dims}} provides convenience methods for accessing
#'   dimensions of pibblefit object
#'   
# Generic functions including \code{\link[=summary.pibblefit]{summary}},
# \code{\link[=print.pibblefit]{print}},
#  \code{\link[=coef.pibblefit]{coef}},
#  \code{\link[=as.list.pibblefit]{as.list}},
#  \code{\link[=predict.pibblefit]{predict}},
#  \code{\link[=model.matrix.pibblefit]{model.matrix}},
#  \code{\link[=name.pibblefit]{name}}, and
#  \code{\link[=sample_prior.pibblefit]{sample_prior}}
#  \code{\link{name_dims}}
# 
# Plotting functions provided by \code{\link[=plot.pibblefit]{plot}}
# and \code{\link[=ppc.pibblefit]{ppc}} (posterior predictive checks)
NULL

#' @rdname orthus_fit
#' @export
#' @references JD Silverman K Roche, ZC Holmes, LA David, S Mukherjee. 
#'   Bayesian Multinomial Logistic Normal Models through Marginally Latent Matrix-T Processes. 
#'   2019, arXiv e-prints, arXiv:1903.11695
orthus <- function(Y=NULL, Z=NULL, X=NULL, upsilon=NULL, Theta=NULL, Gamma=NULL, Xi=NULL,
                   init=NULL, 
                   pars=c("Eta", "Lambda", "Sigma"),
                   ...){
  args <- list(...)
  N <- try_set_dims(c(ncol(Y), ncol(X), args[["N"]]))
  P <- try_set_dims(c(nrow(Z), args[["P"]]))
  D <- try_set_dims(c(nrow(Y), nrow(Theta)+1-P, nrow(Xi)+1-P, ncol(Xi)+1-P, args[["D"]]))
  Q <- try_set_dims(c(nrow(X), ncol(Theta), nrow(Gamma), ncol(Gamma), args[["Q"]]))
  if (any(c(N, D, Q, P) <=0)) stop("N, D, and Q must all be greater than 0 (D must be greater than 1)")
  if (D <= 1) stop("D must be greater than 1")

  ## construct default values ##
  # for priors
  if (is.null(upsilon)) upsilon <- D+P+3  # default is minimal information 
  # but with defined mean
  if (is.null(Theta)) Theta <- matrix(0, D-1+P, Q) # default is mean zero
  if (is.null(Gamma)) Gamma <- diag(Q) # default is iid
  if (is.null(Xi)) {
    ## default is iid on base scale for multinomial parameters and independent for Z dims
    Xi <- diag(D-1+P)
    Xi[1:(D-1), 1:(D-1)] <- matrix(0.5, D-1, D-1)
    diag(Xi) <- 1
    Xi <- Xi * (upsilon-D-P)  # make inverse wishart mean Xi as in previous lines 
  }
  
  # check dimensions
  check_dims(upsilon, 1, "upsilon")
  check_dims(Theta, c(D-1+P, Q), "Theta")
  check_dims(Gamma, c(Q, Q), "Gamma")
  check_dims(Xi, c(D-1+P, D-1+P), "Xi")
  
  # set number of iterations 
  n_samples <- args_null("n_samples", args, 2000)
  use_names <- args_null("use_names", args, TRUE)
  
  # This is the signal to sample the prior only
  if (!(is.null(Y)==is.null(Z))) stop("Y and Z must either both be present or absent")
  if (is.null(Y)){
    if (("Eta" %in% pars) & (is.null(X))) stop("X must be given if Eta is to be sampled")
    stop("sorry but sample_prior.orthusfit is not yet implemented")
    # create orthusfit object and pass to sample_prior then return
    out <- orthusfit(N=N, D=D, Q=Q, P=P, coord_system="alr", alr_base=D, 
                     upsilon=upsilon, Theta=Theta, 
                     Gamma=Gamma, Xi=Xi, 
                     # names_categories=rownames(Y), # these won't be present... 
                     # names_samples=colnames(Y), 
                     # names_covariates=colnames(X), 
                     X=X)
    out <- sample_prior(out, n_samples=n_samples, pars=pars, use_names=use_names)
    return(out)
  } else {
    if (is.null(X)) stop("X must be given to fit model")
    if(is.null(init)) init <- random_pibble_init(Y)   # initialize init should
                                                      # work for orthus too
  }
  
  
  
  # for optimization and laplace approximation
  calcGradHess <- args_null("calcGradHess", args, TRUE)
  b1 <- args_null("b1", args, 0.9)
  b2 <- args_null("b2", args, 0.99)
  step_size <- args_null("step_size", args, 0.003)
  epsilon <- args_null("epsilon", args, 10e-7)
  eps_f <- args_null("eps_f", args, 1e-10)
  eps_g <- args_null("eps_g", args, 1e-4)
  max_iter <- args_null("max_iter", args, 10000)
  verbose <- args_null("verbose", args, FALSE)
  verbose_rate <- args_null("verbose_rate", args, 10)
  decomp_method <- args_null("decomp_method", args, "cholesky")
  eigvalthresh <- args_null("eigvalthresh", args, 0)
  jitter <- args_null("jitter", args, 0)
  multDirichletBoot <- args_null("multDirichletBoot", args, -1.0)
  optim_method <- args_null("optim_method", args, "lbfgs")
  useSylv <- args_null("useSylv", args, TRUE)
  ncores <- args_null("ncores", args, -1)
  seed <- args_null("seed", args, sample(1:2^15, 1))
  
  
  ## precomputation ## 
  # The following is a trick that allows pibble to be used to fit orthus models
  # it relies on collapsing the orthus model to a pibble model using the 
  # conditional form of the matrix-t distribution. 
  A <- diag(N) + t(X) %*% Gamma %*% X
  K <- Xi
  AInv <- chol2inv(chol(A))
  one <- 1:(D-1)
  two <- D:(D-1+P)
  B <- Theta%*%X
  E2 <- Z - B[two,]
  K22Inv <- chol2inv(chol(K[two,two]))
  upsilon.star <- upsilon+P
  K.star <- K[one,one] - K[one,two] %*% K22Inv %*% K[two,one]
  B.star <- B[one,] + K[one,two] %*% K22Inv %*% E2
  A.star <- A%*%(diag(N) + AInv %*% t(E2) %*% K22Inv %*% E2)
  
  # Free up some memory
  rm(A, K, AInv, E2, K22Inv)
  
  if (verbose) cat("Inverting (star) Priors\n")
  K.starInv <- chol2inv(chol(K.star))
  A.starInv <- chol2inv(chol(A.star))
  if (verbose) cat("Starting Optimization\n")
  ## fit collapsed model ##
  fitc <- optimPibbleCollapsed(Y, upsilon.star, B.star, K.starInv, A.starInv, init, n_samples, 
                               calcGradHess, b1, b2, step_size, epsilon, eps_f, 
                               eps_g, max_iter, verbose, verbose_rate, 
                               decomp_method, optim_method, eigvalthresh, 
                               jitter, multDirichletBoot, 
                               useSylv, ncores, seed)
  timerc <- parse_timer_seconds(fitc$Timer)
  
  
  
  # if n_samples=0 or if hessian fails, then use MAP eta estimate for 
  # uncollapsing and unless otherwise specified against, use only the 
  # posterior mean for Lambda and Sigma 
  if (is.null(fitc$Samples)) {
    fitc$Samples <- add_array_dim(fitc$Pars, 3)
    ret_mean <- args_null("ret_mean", args, TRUE)
    if (ret_mean && n_samples>0){
      warning("Laplace Approximation Failed, using MAP estimate of eta", 
              " to obtain Posterior mean of Lambda and Sigma", 
              " (i.e., not sampling from posterior distribution of Lambda or Sigma)")
    }
    if (!ret_mean && n_samples > 0){
      warning("Laplace Approximation Failed, using MAP estimate of eta", 
              "but ret_mean was manually specified as FALSE so sampling", 
              "from posterior of Lambda and Sigma rather than using posterior mean")
    }
  } else {
    ret_mean <- args_null("ret_mean", args, FALSE)
  }
  
  seed <- seed + sample(1:2^15, 1)
  ## uncollapse collapsed model ##
  samples <- array(0, dim=c(D-1+P, N, dim(fitc$Samples)[3]))
  samples[one,,] <- fitc$Samples
  for (i in 1:dim(fitc$Samples)[3]) samples[two,,i] <- Z
  
  fitu <- uncollapsePibble(samples, X, Theta, Gamma, Xi, upsilon, 
                           ret_mean=ret_mean, ncores=ncores, seed=seed)
  timeru <- parse_timer_seconds(fitu$Timer)
  
  timer <- c(timerc, timeru)
  timer <- timer[which(names(timer)!="Overall")]
  timer <- c(timer, 
             "Overall" = unname(timerc["Overall"]) +  unname(timeru["Overall"]), 
             "Uncollapse_Overall" = timeru["Overall"])
  
  
  # Marginal Likelihood Computation
  # not yet implemented for orthus - below is code for pibble
  #d <- D^2 + N*D + D*Q
  #logMarginalLikelihood <- fitc$LogLik+d/2*log(2*pi)+.5*fitc$logInvNegHessDet-d/2*log(N)
  
  
  
  ## pretty output ##
  out <- list()
  if ("Eta" %in% pars){
    out[["Eta"]] <- fitc$Samples
  }
  if ("Lambda" %in% pars){
    out[["Lambda"]] <- fitu$Lambda
  }
  if ("Sigma" %in% pars){
    out[["Sigma"]] <- fitu$Sigma
  }
  
  # By default just returns all other parameters
  out$N <- N
  out$Q <- Q
  out$P <- P
  out$D <- D
  out$Y <- Y
  out$Z <- Z
  out$upsilon <- upsilon
  out$Theta <- Theta
  out$X <- X
  out$Xi <- Xi
  out$Gamma <- Gamma
  out$init <- init
  out$iter <- dim(fitc$Samples)[3]
  # for other methods
  out$names_categories <- rownames(Y)
  out$names_samples <- colnames(Y)
  out$names_covariates <- rownames(X)
  out$names_Zdimensions <- rownames(Z)
  out$coord_system <- "alr"
  out$alr_base <- D
  out$summary <- NULL
  out$Timer <- timer
  #out$logMarginalLikelihood <- logMarginalLikelihood
  attr(out, "class") <- c("orthusfit")
  # add names if present 
  #if (use_names) out <- name(out)
  verify(out) # verify the pibblefit object
  return(out)
}