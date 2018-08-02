% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics_s3.R, R/mongrelfit_methods.R
\name{name_dims}
\alias{name_dims}
\alias{names_covariates}
\alias{names_samples}
\alias{names_categories}
\alias{names_covariates<-}
\alias{names_samples<-}
\alias{names_categories<-}
\alias{names_covariates.mongrelfit}
\alias{names_samples.mongrelfit}
\alias{names_categories.mongrelfit}
\alias{names_covariates<-.mongrelfit}
\alias{names_samples<-.mongrelfit}
\alias{names_categories<-.mongrelfit}
\title{Generic method for getting and setting dimension names of fit object}
\usage{
names_covariates(m)

names_samples(m)

names_categories(m)

names_covariates(m) <- value

names_samples(m) <- value

names_categories(m) <- value

\method{names_covariates}{mongrelfit}(m)

\method{names_samples}{mongrelfit}(m)

\method{names_categories}{mongrelfit}(m)

\method{names_covariates}{mongrelfit}(m) <- value

\method{names_samples}{mongrelfit}(m) <- value

\method{names_categories}{mongrelfit}(m) <- value
}
\arguments{
\item{m}{object}

\item{names}{character vector (or NULL)}
}
\description{
Generic method for getting and setting dimension names of fit object
}