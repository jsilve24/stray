#' #' Just to tie into tidybayes
#' predicted_samples.mongrelfit <- function(...){
#'   predict(...)
#' }

#' @import tidybayes
#' @import ggplot2
plot.mongrelfit <- function(m, par="Lambda", focus.cov=NULL, focus.coord=NULL, 
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
  return(p)
}


#' currently iter must be <= m$iter
ppc <- function(m, iter=50){
  
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