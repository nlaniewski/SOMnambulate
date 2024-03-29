% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/debarcode_funcs.R
\name{get.valley}
\alias{get.valley}
\title{Find a density distribution 'valley'}
\usage{
get.valley(x, quantile.trim = T, no.valley.override = NULL)
}
\arguments{
\item{x}{Numeric vector}

\item{quantile.trim}{Logical. By default, \code{TRUE}; trims 'extreme' values.}

\item{no.valley.override}{Numeric. By default, \code{NULL}; If defined and no valley is detected, this value will be used in place.}
}
\value{
Numeric representing the 'valley' value in a density distribution.
}
\description{
Using the second derivative of the differences in a density distribution's y values, both peaks and valleys can be detected; if a single valley is detected, its x value is returned. A standard deviation adjustment is included to 'adjust/shift' the valley -- used during 'debarcoding' to fine-tune assignment results.
}
