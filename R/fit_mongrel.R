#' Interface to fit mongrel models 
#' 
#' This function is largely a more user friendly wrapper around 
#' \code{\link{optimMongrelCollapsed}} and 
#' \code{\link{uncollapseMongrelCollapsed}}. 
#' See details for model specification. 
#'  Notation: \code{N} is number of samples,
#'  \code{D} is number of multinomial categories, \code{Q} is number
#'  of covariates, \code{iter} is the number of samples of \code{eta} (e.g.,
#'  the parameter \code{n_samples} in the function 
#'  \code{\link{optimMongrelCollapsed}})
#' @param Y D x N matrix of counts
#' @param X Q x N matrix of covariates (design matrix)
#' @param upsilon dof for inverse wishart prior (numeric must be > D) 
#'   (default: D+3)
#' @param Theta (D-1) x Q matrix of prior mean for regression parameters
#'   (default: matrix(0, D-1, Q))
#' @param Gamma QxQ prior covariance matrix 
#'   (default: diag(Q))
#' @param Xi (D-1)x(D-1) prior covariance matrix
#'   (default: ALR transform of diag(1)*(upsilon-D-2)/2 - this is 
#'   essentially iid on "base scale" using Aitchison terminology)
#' @param init (D-1) x Q initialization for regression parameters for 
#'   optimization
#' @param pars character vector of posterior parameters to return
#' @param ... argumenst passed to \code{\link{optimMongrelCollapsed}} and 
#'   \code{\link{uncollapseMongrelCollapsed}}
#' 
#' @details the full model is given by:
#'    \deqn{Y_j ~ Multinomial(Pi_j)}
#'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
#'    \deqn{Eta ~ MN_{D-1 x N}(Lambda*X, Sigma, I_N)}
#'    \deqn{Lambda ~ MN_{D-1 x Q}(Theta, Sigma, Gamma)}
#'    \deqn{Sigma ~ InvWish(upsilon, Xi)}
#'  Where A = (I_N + X * Gamma * X')^{-1}, K^{-1} = Xi is a (D-1)x(D-1) 
#'  covariance matrix, Gamma is a Q x Q covariance matrix, and Phi^{-1} is 
#'  ALRInv_D transform. 
#'  
#'  Default behavior is to use MAP estimate for uncollaping collapsed mongrel 
#'  model if laplace approximation is not preformed. 
#' @return an object of class mongrelfit
#' @md
#' @examples 
#' sim <- mongrel_sim()
#' attach(sim)
#' fit <- fit_mongrel(Y, X)
#' @seealso \code{\link{mongrel_transforms}} provide convenience methods for 
#'  transforming the representation of mongrelfit objects (e.g., conversion to 
#'  proportions, alr, clr, or ilr coordinates.)
#'  
#' \code{\link{access_dims}} provides convenience methods for accessing
#'   dimensions of mongrelfit object
#'   
#' Generic functions including \code{\link[=summary.mongrelfit]{summary}},  
#' \code{\link[=print.mongrelfit]{print}}, 
#'  \code{\link[=coefficients.mongrelfit]{coeficients}},  
#'  \code{\link[=as.list.mongrelfit]{as.list}}, and  
#'  \code{\link[=predict.mongrelfit]{predict}}
#' 
#' Plotting functions provided by \code{\link{plot}} and \code{\link{ppc}} 
#' (posterior predictive checks)
fit_mongrel <- function(Y, X, upsilon=NULL, Theta=NULL, Gamma=NULL, Xi=NULL,
                    init=random_mongrel_init(Y), 
                    pars=c("Eta", "Lambda", "Sigma"),
                    ...){
  args <- list(...)
  
  N <- ncol(Y)
  D <- nrow(Y)
  Q <- nrow(X)
  if (ncol(X) != N) stop("X must have same number of columns as Y")

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
    Xi <- Xi*(upsilon-D-2) # make inverse wishart mean Xi as in previous lines 
  }
  
  # check dimensions
  check_dims(upsilon, 1, "upsilon")
  check_dims(Theta, c(D-1, Q), "Theta")
  check_dims(Gamma, c(Q, Q), "Gamma")
  check_dims(Xi, c(D-1, D-1), "Xi")

  # inline to make default checking of extra arguments easier
  args_null <- function(par, argl, default){
    if (is.null(argl[[par]])) return(default)
    return(argl[[par]])
  }
  
  # for optimization and laplace approximation
  n_samples <- args_null("n_samples", args, 2000)
  calcGradHess <- args_null("calcGradHess", args, TRUE)
  b1 <- args_null("b1", args, 0.9)
  b2 <- args_null("b2", args, 0.99)
  step_size <- args_null("step_size", args, 0.003)
  epsilon <- args_null("epsilon", args, 10e-7)
  eps_f <- args_null("eps_f", args, 1e-8)
  eps_g <- args_null("eps_g", args, 1e-5)
  max_iter <- args_null("max_iter", args, 10000)
  verbose <- args_null("verbose", args, FALSE)
  verbose_rate <- args_null("verbose_rate", args, 10)
  decomp_method <- args_null("decomp_method", args, "eigen")
  eigvalthresh <- args_null("eigvalthresh", args, 0)
  
  ## precomputation ## 
  # consider moving these to optim function in the future for numerical
  # stability
  K <- solve(Xi)
  A <- solve(diag(N) + t(X) %*% Gamma %*% X)

  ## fit collapsed model ##
  fitc <- optimMongrelCollapsed(Y, upsilon, Theta%*%X, K, A, init, n_samples, 
                                calcGradHess, b1, b2, step_size, epsilon, eps_f, 
                                eps_g, max_iter, verbose, verbose_rate, 
                                decomp_method, eigvalthresh, no_error=TRUE)
  # if n_samples=0 or if hessian fails, then use MAP eta estimate for 
  # uncollapsing and unless otherwise specified against, use only the 
  # posterior mean for Lambda and Sigma 
  if (is.null(fitc$Samples)) {
    fitc$Samples <- fitc$pars
    ret_mean <- args_null("ret_mean", args, TRUE)
    if (ret_mean && n_samples>0){
      warning("Laplace Approximation Failed, using MAP estimate of eta", 
              "to obtain Posterior mean of Lambda and Sigma", 
              "i.e., not sampling from posterior distribution of Lambda or Sigma")
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
  fitu <- uncollapseMongrelCollapsed(fitc$Samples, X, Theta, Gamma, Xi, upsilon, 
                                     ret_mean=ret_mean)
  
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
  
  attr(out, "class") <- c("mongrelfit")
  # add names if present 
  out <- apply_names_mongrel(out)
  
  return(out)
}