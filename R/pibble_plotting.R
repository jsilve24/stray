#' Plot Summaries of Posterior Distribution of pibblefit Parameters
#' 
#' @param x an object of class pibblefit
#' @param ... other arguments passed to plot.pibblefit (see details)
#' @return ggplot object
#' @import tidybayes
#' @import ggplot2
#' @importFrom dplyr filter
#' @details Other arguments:
#' \itemize{
#' \item `par` parameter to plot (options: Lambda, Eta, and Sigma) 
#'   (default="Lambda")
#' \item `focus.cov` vector of covariates to include in plot (plots all if NULL)
#' \item `focus.coord` vector of coordinates to include in plot (plots all if NULL)
#' \item `focus.sample` vector of samples to include in plot (plots all if NULL)
#' \item `use_names` if TRUE, uses dimension names found in data as plot labels
#'   rather than using dimension integer indices. 
#' }
#' @export
#' @examples
#' sim <- pibble_sim(N=10, D=4, Q=3)
#' fit <- pibble(sim$Y, sim$X)
#' plot(fit, par="Lambda")
#' plot(fit, par="Sigma")
plot.pibblefit <- function(x, ...){
 args <- list(...)
 par <- args_null("par", args, "Lambda")
 focus.cov <- args_null("focus.cov", args, NULL)
 focus.coord <- args_null("focus.coord", args, NULL)
 focus.sample<- args_null("focus.sample", args, NULL)
 use_names <- args_null("use_names", args, TRUE)

 if (is.null(x[[par]])) stop("pibblefit object does not contain samples for specified parameter")
 if (par %in% c("Lambda", "Eta")) {
   return(plot_mf_lambdaeta(x, par, focus.cov, focus.coord, 
                            focus.sample, use_names))
 } else if (par=="Sigma"){
   return(plot_mf_sigma(x, focus.coord, use_names))
 } else {
   stop("only parameters Lambda, Sigma, and Eta are recognized by plot.pibblefit")
 }
}

plot_mf_lambdaeta <- function(m, par, focus.cov=NULL, focus.coord=NULL, 
                              focus.sample=NULL, use_names=TRUE){
  data <- summary.pibblefit(m, pars=par, use_names=use_names, 
                             as_factor=TRUE, gather_prob=TRUE)
  data <- data[[par]]
    
  # some code to handle numeric focuses
  
  # Focus 
  if (!is.null(focus.cov)) data <- filter(data, .data$covariate %in% focus.cov)
  if (!is.null(focus.coord)) data <- filter(data, .data$coord %in% focus.coord)
  if (!is.null(focus.sample)) data <- filter(data, .data$sample %in% focus.sample)
  
  
  if (par=="Lambda"){
    p <- data %>% 
      ggplot(aes(x=.data$val, y=.data$coord)) +
      geom_intervalh() +
      geom_point() + 
      facet_grid(~.data$covariate)
  }
  if (par=="Eta"){
    p <- data %>% 
      ggplot(aes(x=.data$val, y=.data$coord)) +
      geom_intervalh() +
      geom_point() +
      facet_grid(~.data$sample)
  }
  p <- p+
    theme_minimal() +
    scale_color_brewer() +
    guides(color=guide_legend(title="Credible\nInterval")) +
    theme(axis.title.y=element_blank())
  # Set axis labels
  if (m$coord_system %in% c("clr", "ilr", "alr")){
    p <- p + xlab("Log-Ratio Value")
  } else if (m$coord_system == "proportions"){
    p <- p + xlab("Proportions")
  }
  
  return(p)
}

plot_mf_sigma <- function(m, focus.coord=NULL, use_names=TRUE){
  data <- pibble_tidy_samples(m, use_names, as_factor=TRUE)
  if (!is.null(focus.coord)) data <- filter(data, 
                                            .data$coord %in% focus.coord, 
                                            .data$coord2 %in% focus.coord)
  data <- dplyr::filter(data, .data$Parameter=="Sigma")
  data <-  group_by(data, .data$coord, .data$coord2) %>% 
    mutate(medval = median(.data$val)) %>% 
    ungroup() 
  p <- data %>% 
    ggplot(aes(x = .data$val)) +
    geom_rect(data = filter(data, .data$iter==1), 
              aes(fill=.data$medval), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  
  if (m$iter > 1) p <- p+ geom_density(fill="lightgrey") 
  p <- p +
    facet_grid(.data$coord~.data$coord2) +
    theme_minimal() +
    theme(strip.text.y = element_text(angle=0)) +
    scale_fill_distiller(name="Median Value", palette = "Blues") +
    xlab("Value") +
    ylab("Density")
  return(p)
}



#' Visualization of Posterior Predictive Check of fit model
#' 
#' @param m an object of class pibblefit
#' @param ... other options passed to ppc (see details)
#' @return ggplot object
#' @import ggplot2
#' @importFrom driver gather_array
#' @importFrom stats quantile
#' @export
#' @details ppc.pibblefit accepts the following additional arguments:
#' \itemize{
#'   \item "type" type of plot (options "lines",  "points", "bounds")
#'   \item "iter" number of samples from posterior predictive distribution to plot
#'   (currently must be <= m$iter) if type=="lines" default is 50, if type=="ribbon"
#'   default is to use all available iterations. 
#'   \item "from_scratch" should predictions of Y come from fitted Eta or from 
#'   predictions of Eta from posterior of Lambda? (default: false)
#' }
#' @examples 
#' \dontrun{
#' fit <- pibble(Y, X)
#' ppc(fit)
#' }
ppc.pibblefit <- function(m, ...){
  args <- list(...)
  type <- args_null("type", args, "bounds")
  iter <-  args_null("iter", args, NULL)
  from_scratch <- args_null("from_scratch", args, FALSE)
  
  msg <- paste("No observed count data (Y) to check against", 
               "perhaps you are looking for the function `predict`?")
  if (is.null(m$Y)) stop(msg)
  
  if (is.null(iter)) {
    if (type =="lines") iter <- min(50, niter(m))
    if (type %in% c("bounds", "points")) iter <- niter(m)
  } 
  if (iter > m$iter) {
    stop("iter param larger than number of stored posterior samples")
  }
  
  if (!is.null(m$Y)) o <- order(m$Y, decreasing=TRUE) else o <- 1:(m$N*m$D)
  
  pp <- predict(m, response="Y", from_scratch=from_scratch)
  if (iter < niter(m)) pp <- pp[,,sample(1:m$iter, iter), drop=F]
  pp <- matrix(pp, m$D*m$N, iter) 
  pp <- pp[o,]
  
  tr <- data.frame(dim_1 = 1:(m$N*m$D), 
                   dim_2 = NA, 
                   val = c(m$Y)[o])
  
  if (type %in% c("bounds", "points")){
    pp <- apply(pp, 1, function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.975)))
    rownames(pp) <- c("p2.5", "p50", "p97.5")
    pp <- as.data.frame(t(pp))
    pp$dim_1 <- 1:nrow(pp)
    p <- ggplot(pp, aes(x=.data$dim_1, y=.data$p50)) 
    if (type == "bounds") {
      p <- p + geom_linerange(aes(ymin=.data$p2.5, ymax=.data$p97.5), color="black", alpha=0.4)
    }
     p <- p + 
       geom_point(color="lightgrey", size=1) +
       geom_line(data=tr, aes(y=.data$val), color="green", alpha=0.6)
  } else if (type == "lines"){
    pp <- driver::gather_array(pp, .data$val) 
    p <- ggplot(pp, aes(x = .data$dim_1, y = .data$val)) +
      geom_line(aes(group=.data$dim_2), color="black", alpha=0.4) +
      geom_line(data=tr, color="green", alpha=0.6)
  }
  p <- p + theme_minimal() +xlab("") +ylab("Counts")
  return(p)
}
