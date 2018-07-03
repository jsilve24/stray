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
  if (use_names){
    if (!is.null(m$names_samples)) {
      l <- dplyr::mutate(l, sample = m$names_samples[sample])
    }
    if (!is.null(m$names_covariates)){
      l <- dplyr::mutate(l, covariate = m$names_covariates[covariate])
    }
    if (!is.null(m$names_categories)){
      if (m$coord_system == "proportions"){
        l <- dplyr::mutate(l, coord = paste0("prop_", m$names_categories[coord]))
      } 
      if (m$coord_system == "clr"){
        l <- dplyr::mutate(l, coord = paste0("clr_", m$names_categories[coord]), 
                    coord2 = paste0("clr_", m$names_categories[coord2]))
      } 
    }
  }
  return(l)
}

#' @param ... other expressions to pass to summarise (using name 'val' unquoted is 
#'   probably what you want)
#' @import dplyr
#' @importFrom driver summarise_posterior
#' @examples 
#' ... what is fit ... 
#' summary(fit, pars="Eta", median = median(val))
summary.mongrelfit <- function(m, pars=NULL, use_names=TRUE, ...){
  if (is.null(pars)) pars <- c("Eta", "Lambda", "Sigma")
  mtidy <- dplyr::filter(mongrel_tidy_samples(m, use_names), Parameter %in% pars)
  if (m$coord_system != "proportions") {
    mtidy %>%
      group_by(Parameter, coord, coord2, sample, covariate) %>% 
      summarise_posterior(val, ...) %>%
      ungroup() %>%
      return()
  } else {
    mtidy %>%
      # coord2 does not exist if in proportions 
      group_by(Parameter, coord, sample, covariate) %>%
      summarise_posterior(val, ...) %>%
      ungroup() %>%
      return()
  }
}

# print should give dimensions and coordinate system thing is in.  
# plot




