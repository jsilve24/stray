#' Interface to fit pibble models 
#' 
#' This function is largely a more user friendly wrapper around 
#' \code{\link{optimPibbleCollapsed}} and 
#' \code{\link{uncollapsePibble}}. 
#' See details for model specification. 
#'  Notation: \code{N} is number of samples,
#'  \code{D} is number of multinomial categories, \code{Q} is number
#'  of covariates, \code{iter} is the number of samples of \code{eta} (e.g.,
#'  the parameter \code{n_samples} in the function 
#'  \code{\link{optimPibbleCollapsed}})
#' @param Y D x N matrix of counts (if NULL uses priors only)
#' @param X Q x N matrix of covariates (design matrix) (if NULL uses priors only, must
#' be present to sample Eta)
#' @param upsilon dof for inverse wishart prior (numeric must be > D) 
#'   (default: D+3)
#' @param Theta (D-1) x Q matrix of prior mean for regression parameters
#'   (default: matrix(0, D-1, Q))
#' @param Gamma QxQ prior covariance matrix 
#'   (default: diag(Q))
#' @param Xi (D-1)x(D-1) prior covariance matrix
#'   (default: ALR transform of diag(1)*(upsilon-D)/2 - this is 
#'   essentially iid on "base scale" using Aitchison terminology)
#' @param init (D-1) x Q initialization for Eta for optimization
#' @param pars character vector of posterior parameters to return
#' @param m object of class pibblefit 
#' @param ... arguments passed to \code{\link{optimPibbleCollapsed}} and 
#'   \code{\link{uncollapsePibble}}
#' 
#' @details the full model is given by:
#'    \deqn{Y_j \sim Multinomial(Pi_j)} 
#'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
#'    \deqn{Eta \sim MN_{D-1 x N}(Lambda*X, Sigma, I_N)}
#'    \deqn{Lambda \sim MN_{D-1 x Q}(Theta, Sigma, Gamma)}
#'    \deqn{Sigma \sim InvWish(upsilon, Xi)}
#'  Where Gamma is a Q x Q covariance matrix, and Phi^{-1} is 
#'  ALRInv_D transform. 
#'  
#'  Default behavior is to use MAP estimate for uncollaping the LTP 
#'  model if laplace approximation is not preformed. 
#' @return an object of class pibblefit
#' @md
#' @name pibble_fit
#' @examples 
#' sim <- pibble_sim()
#' fit <- pibble(sim$Y, sim$X)
#' @seealso \code{\link{stray_transforms}} provide convenience methods for 
#'  transforming the representation of pibblefit objects (e.g., conversion to 
#'  proportions, alr, clr, or ilr coordinates.)
#'  
#' \code{\link{access_dims}} provides convenience methods for accessing
#'   dimensions of pibblefit object
#'   
#' Generic functions including \code{\link[=summary.pibblefit]{summary}},  
#' \code{\link[=print.pibblefit]{print}}, 
#'  \code{\link[=coef.pibblefit]{coef}},  
#'  \code{\link[=as.list.pibblefit]{as.list}},  
#'  \code{\link[=predict.pibblefit]{predict}}, 
#'  \code{\link[=model.matrix.pibblefit]{model.matrix}},
#'  \code{\link[=name.pibblefit]{name}}, and
#'  \code{\link[=sample_prior.pibblefit]{sample_prior}}
#'  \code{\link{name_dims}}
#' 
#' Plotting functions provided by \code{\link[=plot.pibblefit]{plot}} 
#' and \code{\link[=ppc.pibblefit]{ppc}} (posterior predictive checks)
NULL

#' @rdname pibble_fit
#' @export
#' @references JD Silverman K Roche, ZC Holmes, LA David, S Mukherjee. 
#'   Bayesian Multinomial Logistic Normal Models through Marginally Latent Matrix-T Processes. 
#'   2019, arXiv e-prints, arXiv:1903.11695
pibble <- function(Y=NULL, X=NULL, upsilon=NULL, Theta=NULL, Gamma=NULL, Xi=NULL,
                    init=NULL, 
                    pars=c("Eta", "Lambda", "Sigma"),
                    ...){
  args <- list(...)
  N <- try_set_dims(c(ncol(Y), ncol(X), args[["N"]]))
  D <- try_set_dims(c(nrow(Y), nrow(Theta)+1, nrow(Xi)+1, ncol(Xi)+1, args[["D"]]))
  Q <- try_set_dims(c(nrow(X), ncol(Theta), nrow(Gamma), ncol(Gamma), args[["Q"]]))
  if (any(c(N, D, Q) <=0)) stop("N, D, and Q must all be greater than 0 (D must be greater than 1)")
  if (D <= 1) stop("D must be greater than 1")
  
  ## construct default values ##
  # for priors
  if (is.null(upsilon)) upsilon <- D+3  # default is minimal information 
                                        # but with defined mean
  if (is.null(Theta)) Theta <- matrix(0, D-1, Q) # default is mean zero
  if (is.null(Gamma)) Gamma <- diag(Q) # default is iid
  if (is.null(Xi)) {
    # default is iid on base scale
    # G <- cbind(diag(D-1), -1) ## alr log-constrast matrix
    # Xi <- 0.5*G%*%diag(D)%*%t(G) ## default is iid on base scale
    Xi <- matrix(0.5, D-1, D-1) # same as commented out above 2 lines
    diag(Xi) <- 1               # same as commented out above 2 lines
    Xi <- Xi*(upsilon-D) # make inverse wishart mean Xi as in previous lines 
  }
  
  # check dimensions
  check_dims(upsilon, 1, "upsilon")
  check_dims(Theta, c(D-1, Q), "Theta")
  check_dims(Gamma, c(Q, Q), "Gamma")
  check_dims(Xi, c(D-1, D-1), "Xi")
  
  # set number of iterations 
  n_samples <- args_null("n_samples", args, 2000)
  use_names <- args_null("use_names", args, TRUE)
  
  # This is the signal to sample the prior only
  if (is.null(Y)){
    if (("Eta" %in% pars) & (is.null(X))) stop("X must be given if Eta is to be sampled")
    # create pibblefit object and pass to sample_prior then return
    out <- pibblefit(N=N, D=D, Q=Q, coord_system="alr", alr_base=D, 
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
    if(is.null(init)) init <- random_pibble_init(Y)   # initialize init 
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
  

  ## precomputation ## 
  #KInv <- solve(Xi)
  #AInv <- solve(diag(N) + t(X) %*% Gamma %*% X)
  if (verbose) cat("Inverting Priors\n")
  KInv <- chol2inv(chol(Xi))
  AInv <- chol2inv(chol(diag(N) + t(X) %*% Gamma %*% X))
  if (verbose) cat("Starting Optimization\n")
  ## fit collapsed model ##
  fitc <- optimPibbleCollapsed(Y, upsilon, Theta%*%X, KInv, AInv, init, n_samples, 
                                calcGradHess, b1, b2, step_size, epsilon, eps_f, 
                                eps_g, max_iter, verbose, verbose_rate, 
                                decomp_method, optim_method, eigvalthresh, 
                                jitter, multDirichletBoot, 
                                useSylv, ncores)
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
  
  seed <- args_null("seed", args, sample(1:2^15, 1))
  ## uncollapse collapsed model ##
  fitu <- uncollapsePibble(fitc$Samples, X, Theta, Gamma, Xi, upsilon, 
                                     ret_mean=ret_mean, ncores=ncores, seed=seed)
  timeru <- parse_timer_seconds(fitu$Timer)
  
  timer <- c(timerc, timeru)
  timer <- timer[which(names(timer)!="Overall")]
  timer <- c(timer, 
             "Overall" = unname(timerc["Overall"]) +  unname(timeru["Overall"]), 
             "Uncollapse_Overall" = timeru["Overall"])
  
  
  # Marginal Likelihood Computation
  d <- D^2 + N*D + D*Q
  logMarginalLikelihood <- fitc$LogLik+d/2*log(2*pi)+.5*fitc$logInvNegHessDet-d/2*log(N)
  
  
  
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
  out$D <- D
  out$Y <- Y
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
  out$coord_system <- "alr"
  out$alr_base <- D
  out$summary <- NULL
  out$Timer <- timer
  out$logMarginalLikelihood <- logMarginalLikelihood
  attr(out, "class") <- c("pibblefit")
  # add names if present 
  if (use_names) out <- name(out)
  verify(out) # verify the pibblefit object
  return(out)
}

#' @title mongrel
#' @inheritParams pibble
#' @description This function is deprecated, please use \code{pibble} 
#' instead. 
#' @name mongrel-deprecated
mongrel <- function(Y=NULL, X=NULL, upsilon=NULL, Theta=NULL, Gamma=NULL, Xi=NULL,
                    init=NULL, 
                    pars=c("Eta", "Lambda", "Sigma"),
                    ...){
  .Deprecated("pibble")
  pibble(Y, X, upsilon, Theta, Gamma, Xi, init, pars, ...)
}



#' @rdname pibble_fit
#' @export
refit.pibblefit <- function(m, pars=c("Eta", "Lambda", "Sigma"), ...){
  # Store coordinates and tranfsorm to cannonical representation
  l <- store_coord(m)
  m <- to_alr(m, m$D)
  
  # Concatenate parameters to pass to pibble function
  argl <- list(...)
  argl$pars <- pars
  ml <- as.list(m)
  argl <- c(ml, argl)
  
  # Need to handle iter as part of m but no n_samples passed
  # in this situation should pull iter from m and pass as n_samples to pibble 
  if (is.null(argl[["n_samples"]]) & !is.null(m$iter)) argl[["n_samples"]] <- m$iter 
  
  # pass to pibble function
  m <- do.call(pibble, argl)
  
  # Reapply original coordinates
  m <- reapply_coord(m, l)
  verify(m)
  return(m)
}
