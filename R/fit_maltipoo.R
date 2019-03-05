#' Interface to fit maltipoo models
#' 
#' This function is largely a more user friendly wrapper around 
#' \code{\link{optimMaltipooCollapsed}} and 
#' \code{\link{uncollapseMongrelCollapsed}}. 
#' See details for model specification. 
#'  Notation: \code{N} is number of samples,
#'  \code{D} is number of multinomial categories, \code{Q} is number
#'  of covariates, \code{P} is the number of variance components 
#'  \code{iter} is the number of samples of \code{eta} (e.g.,
#'  the parameter \code{n_samples} in the function 
#'  \code{\link{optimMongrelCollapsed}})
#'  
#' @param U a PQ x Q matrix of stacked variance components (each of dimension Q x Q)
#' @param init (D-1) x Q initialization for Eta for optimization
#' @param ellinit P vector initialization values for ell for optimization 
#' @inheritParams mongrel_fit
#'  
#' @details the full model is given by:
#'    \deqn{Y_j ~ Multinomial(Pi_j)}
#'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
#'    \deqn{Eta ~ MN_{D-1 x N}(Lambda*X, Sigma, I_N)}
#'    \deqn{Lambda ~ MN_{D-1 x Q}(Theta, Sigma, Gamma)}
#'    \deqn{Gamma = e^{ell_1} U_1 + ... + e^{ell_P} U_P}
#'    \deqn{Sigma ~ InvWish(upsilon, Xi)}
#'    
#'  Where A = (I_N + X * Gamma * X')^{-1}, K^{-1} = Xi is a (D-1)x(D-1) 
#'  covariance matrix, U_1 is a Q x Q covariance matrix (a variance component), 
#'  e^{ell_i} is a scale for that variance component and Phi^{-1} is 
#'  ALRInv_D transform. 
#'  
#'  Default behavior is to use MAP estimate for uncollaping collapsed mongrel 
#'  model if laplace approximation is not preformed. 
#'  
#'  Parameters ell are treated as fixed and estimated by MAP estimation. 
#'  
#' @name maltipoo_fit
#' @return an object of class maltipoofit
NULL

#' @rdname maltipoo_fit
#' @export
maltipoo <- function(Y=NULL, X=NULL, upsilon=NULL, Theta=NULL, U=NULL, 
                     Xi=NULL, init=NULL, ellinit=NULL, 
                     pars=c("Eta", "Lambda", "Sigma"), 
                     ...){
  args <- list(...)
  
  N <- try_set_dims(c(ncol(Y), ncol(X), args[["N"]]))
  D <- try_set_dims(c(nrow(Y), nrow(Theta)+1, nrow(Xi)+1, ncol(Xi)+1, args[["D"]]))
  Q <- try_set_dims(c(nrow(X), ncol(Theta), ncol(U), args[["Q"]]))
  P <- try_set_dims(c(nrow(U)/Q)) # Maltipoo specific
  if (any(c(N, D, Q) <=0)) stop("N, D, and Q must all be greater than 0 (D must be greater than 1)")
  if (D <= 1) stop("D must be greater than 1")
  
  ## construct default values ##
  # for priors
  if (is.null(upsilon)) upsilon <- D+3  # default is minimal information 
  # but with defined mean
  if (is.null(Theta)) Theta <- matrix(0, D-1, Q) # default is mean zero
  if (is.null(U)) U <- diag(Q) # default is iid
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
  check_dims(U, c(P*Q, Q), "Q")
  check_dims(Xi, c(D-1, D-1), "Xi")
  
  # set number of iterations 
  n_samples <- args_null("n_samples", args, 2000)
  use_names <- args_null("use_names", args, TRUE)
  
  # This is the signal to sample the prior only
  if (is.null(Y)){
    if (("Eta" %in% pars) & (is.null(X))) stop("X must be given if Eta is to be sampled")
    # create matipoofit object and pass to sample_prior then return
    out <- maltipoofit(N=N, D=D, Q=Q, P=P, coord_system="alr", alr_base=D, 
                      upsilon=upsilon, Theta=Theta, Xi=Xi,U=U, 
                      # names_categories=rownames(Y), # these won't be present... 
                      # names_samples=colnames(Y), 
                      # names_covariates=colnames(X), 
                      X=X)
    out <- sample_prior(out, n_samples=n_samples, pars=pars, use_names=use_names)
    return(out)
  } else {
    if (is.null(X)) stop("X must be given to fit model")
    if(is.null(init)) init <- random_mongrel_init(Y)   # initialize init 
    if(is.null(ellinit)) ellinit <- rep(0, P)
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
  calcPartialHess <- args_null("calcPartialHess", args, FALSE)
  
  ## precomputation ## 
  K <- solve(Xi)
  
  ## fit collapsed model ##
  fitc <- optimMaltipooCollapsed(Y, upsilon, Theta, X, K, U, init, ellinit, 
                                n_samples, 
                                calcGradHess, b1, b2, step_size, epsilon, eps_f, 
                                eps_g, max_iter, verbose, verbose_rate, 
                                decomp_method, eigvalthresh, 
                                jitter, calcPartialHess)
  
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
  
  ## uncollapse collapsed model ##
  # Calculate Gamma using MAP estimate for ell
  Gamma <- matrix(0, Q, Q)
  for (i in 1:P){
    Gamma <- Gamma + fitc$VCScale[i]*U[((i-1)*Q+1):(i*Q),]
  }
  seed <- args_null("seed", args, sample(1:2^15, 1))
  
  fitu <- uncollapseMongrelCollapsed(fitc$Samples, X, Theta, Gamma, Xi, upsilon, 
                                     ret_mean=ret_mean, seed=seed)
  
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
  
  # Marginal Likelihood
  d <- D^2 + N*D + D*Q + length(fitc$VCScale)
  logMarginalLikelihood <- fitc$LogLik+d/2*log(2*pi)+.5*fitc$logInvNegHessDet - d/2*log(N)
  
  # By default just returns all other parameters
  out$N <- N
  out$Q <- Q
  out$D <- D
  out$P <- P
  out$Y <- Y
  out$upsilon <- upsilon
  out$Theta <- Theta
  out$X <- X
  out$Xi <- Xi
  out$U <- U
  out$VCScale <- fitc$VCScale
  out$init <- init
  out$ellinit <- ellinit
  out$iter <- dim(fitc$Samples)[3]
  # for other methods
  out$names_categories <- rownames(Y)
  out$names_samples <- colnames(Y)
  out$names_covariates <- rownames(X)
  out$coord_system <- "alr"
  out$alr_base <- D
  out$summary <- NULL
  out$logMarginalLikelihood
  attr(out, "class") <- c("maltipoofit", "mongrelfit")
  # add names if present 
  if (use_names) out <- name(out)
  verify(out) # verify the mongrelfit object
  return(out)
}


# Refit
# Sample prior