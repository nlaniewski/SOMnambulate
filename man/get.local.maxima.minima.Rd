% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/debarcode_funcs.R
\name{get.local.maxima.minima}
\alias{get.local.maxima.minima}
\title{Get local maxima and minima from a density distribution}
\usage{
get.local.maxima.minima(x, quantile.probs = NULL, plot = F)
}
\arguments{
\item{x}{Numeric vector.}

\item{quantile.probs}{Numeric vector of length 2; used to define lower/upper quantile \link[stats:quantile]{probabilities}; trims 'extreme events' before calculating the density distribution.}

\item{plot}{Logical; if \code{TRUE}, a plot will be drawn displaying a \link[stats]{density} distribution with vertical \link[graphics]{abline}s: "purple peaks" and "violet valleys".}
}
\value{
a list; contains numeric values: local maxima \verb{[['lmaxima']]} and local minima \verb{[['lminima']]}.
}
\description{
Get local maxima (peaks) and minima (valleys) from a density distribution
}
\examples{
dt<-SOMnambulate:::prepared.examples('dt')

dt[,CD45_112Cd := asinh(CD45_112Cd/10)]
x<-dt[['CD45_112Cd']]

#maxima and minima values returned for extreme tails
SOMnambulate:::get.local.maxima.minima(x,plot=TRUE)

#quantile 'trim'
SOMnambulate:::get.local.maxima.minima(x,quantile.probs=c(0,0.999),plot=TRUE)
}
