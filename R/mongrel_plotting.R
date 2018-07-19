#' Plot Summaries of Posterior Distribution of Mongrelfit Parameters
#' 
#' @param m an object of class mongrelfit
#' @param par parameter to plot (options: Lambda, Eta, and Sigma) 
#'   (default="Lambda")
#' @param focus.cov vector of covariates to include in plot (plots all if NULL)
#' @param focus.coord vector of coordinates to include in plot (plots all if NULL)
#' @param focus.sample vector of samples to include in plot (plots all if NULL)
#' @param use_names if TRUE, uses dimension names found in data as plot labels
#'   rather than using dimension integer indicies. 
#' 
#' @return ggplot object
#' @import tidybayes
#' @import ggplot2
#' @importFrom dplyr filter
#' @export
#' @examples
#' \dontrun{
#' fit <- mongrel(Y, X)
#' plot(fit, par="Lambda")
#' plot(fit, par="Sigma")
#' plot(fit, par="Sigma", focus.coord=c("s1", "s2", "s3"))
#' }
plot.mongrelfit <- function(m, par="Lambda", focus.cov=NULL, focus.coord=NULL, 
                            focus.sample=NULL, use_names=TRUE){
 if (is.null(m[[par]])) stop("mongrelfit object does not contain samples for specified parameter")
 if (par %in% c("Lambda", "Eta")) {
   return(plot_mf_lambdaeta(m, par, focus.cov, focus.coord, 
                            focus.sample, use_names))
 } else if (par=="Sigma"){
   return(plot_mf_sigma(m, focus.coord, use_names))
 } else {
   stop("only parameters Lambda, Sigma, and Eta are recognized by plot.mongrelfit")
 }
}

plot_mf_lambdaeta <- function(m, par, focus.cov=NULL, focus.coord=NULL, 
                              focus.sample=NULL, use_names=TRUE){
  data <- summary.mongrelfit(m, pars=par, use_names=use_names, gather_prob=TRUE)
  data <- data[[par]]
  
  
  # some code to handle numeric focuses
  
  # Focus 
  if (!is.null(focus.cov)) data <- filter(data, covariate %in% focus.cov)
  if (!is.null(focus.coord)) data <- filter(data, coord %in% focus.coord)
  if (!is.null(focus.sample)) data <- filter(data, sample %in% focus.sample)
  
  if (par=="Lambda"){
    p <- data %>% 
      ggplot(aes(x=val, y=coord)) +
      geom_intervalh() +
      geom_point() + 
      facet_grid(~covariate)
  }
  if (par=="Eta"){
    p <- data %>% 
      ggplot(aes(x=val, y=coord)) +
      geom_intervalh() +
      geom_point() +
      facet_grid(~sample)
  }
  p <- p+
    theme_minimal() +
    scale_color_brewer()
    guides(color=guide_legend(title="Credible Interval"))
  return(p)
}

plot_mf_sigma <- function(m, focus.coord=NULL, use_names=TRUE){
  data <- mongrel_tidy_samples(m, use_names)
  if (!is.null(focus.coord)) data <- filter(data, 
                                            coord %in% focus.coord, 
                                            coord2 %in% focus.coord)
  data <-  group_by(data, coord, coord2) %>% 
    mutate(medval = median(val)) %>% 
    ungroup() 
  p <- data %>% 
    ggplot(aes(x = val)) +
    geom_rect(data = filter(data, iter==1), 
              aes(fill=medval), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_density(fill="lightgrey") +
    facet_grid(coord~coord2) +
    theme_minimal() +
    theme(strip.text.y = element_text(angle=0)) +
    scale_fill_distiller(name="Median Value", palette = "Blues") +
    xlab("Value") +
    ylab("Density")
  return(p)
}



#' Visualization of Posterior Predictive Check of fit model
#' 
#' @param m an object of class mongrelfit
#' @param iter number of samples from posterior predictive distribution to plot
#'   (currently must be <= m$iter)
#' @return ggplot object
#' @import ggplot2
#' @importFrom driver gather_array
#' @export
#' @examples 
#' \dontrun{
#' fit <- mongrel(Y, X)
#' ppc(fit)
#' }
ppc.mongrelfit <- function(m, iter=50){
  
  pp <- predict(m, response="Y")
  pp <- pp[,,sample(1:m$iter, iter)]
  pp <- matrix(pp, m$D*m$N, iter) 
  pp <- gather_array(pp, val) 
  tr <- data.frame(dim_1 = 1:(m$N*m$D), 
                   dim_2 = NA, 
                   val = c(m$Y))
  p <- ggplot(pp, aes(x = dim_1, y = val)) +
    geom_line(aes(group=dim_2), color="black", alpha=0.4) +
    geom_line(data=tr, color="green", alpha=0.6)

  p <- p + theme_minimal()
  
  return(p)
}