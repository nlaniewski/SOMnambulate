% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsom_funcs.R
\name{fsom.merge.codes}
\alias{fsom.merge.codes}
\title{Merge FlowSOM codes into a new cluster}
\usage{
fsom.merge.codes(fsom, codes.to.merge)
}
\arguments{
\item{fsom}{FlowSOM object, as returned from \code{FlowSOM::SOM} and previously clustered using \code{fsom.codes.to.clusters}.}

\item{codes.to.merge}{Numeric vector; individual codes to be manually merged into a new cluster.}
}
\value{
FlowSOM object with an updated cluster factor (\code{fsom$cluster$fac})
}
\description{
Merge FlowSOM codes into a new cluster
}
