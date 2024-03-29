% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotfc}
\alias{plotfc}
\title{Title forecast plotting}
\usage{
plotfc(
  fcout,
  ydata,
  dateseq,
  vnames,
  fulldates,
  frequency = 12,
  ymat = c(4.35, 4.7, 4.55, 4.65, 8.2, 8.5, 7, 7.4, 7.2, 7.6, -0.02, 0.06, 5.6, 6.4,
    -0.01, 0.05, 0, 0.12, -0.01, 0.05),
  start_date = NULL,
  filename = "fcout",
  cushion = 0
)
}
\arguments{
\item{fcout}{forecasting results}

\item{ydata}{Original variable dataframe}

\item{dateseq}{Vector of points from which to draw forecasting lines}

\item{vnames}{Names of variables to use for each plot}

\item{fulldates}{full range of dates associated with time series}

\item{frequency}{how often time points are sampled per year (e.g. 12 for
monthly data, 4 for quarterly data)}

\item{cushion}{}
}
\description{
Title forecast plotting
}
