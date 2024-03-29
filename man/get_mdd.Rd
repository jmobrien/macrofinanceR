% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bpss_tools.R
\name{get_mdd}
\alias{get_mdd}
\title{Title Small wrapper function for calculated mdd with MHM methods}
\usage{
get_mdd(
  x,
  lh,
  trunc = 0.95,
  delta = NULL,
  alpha = 3,
  beta = 3,
  efac = NULL,
  covsub = FALSE
)
}
\arguments{
\item{x}{draws in a matrix}

\item{lh}{vector of -log posterior densities (conditional on certain params)}

\item{trunc}{extra piece numerator density for MHM}

\item{delta}{delta_it in an array, if appropriate}

\item{alpha, beta}{parameters for the distribution of the delta (maybe you should change this)}

\item{efac}{extra piece numerator density for MHM}

\item{covsub}{whether to calc covariance matrix with all xout, or 5000 evenly spaced draws}
}
\description{
Title Small wrapper function for calculated mdd with MHM methods
}
