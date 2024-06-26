% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/showmeSOMething.R
\name{showmeSOMething}
\alias{showmeSOMething}
\title{Interactive exploration of high-dimensional clustering (FlowSOM) results}
\usage{
showmeSOMething(
  dat,
  cluster.col = NULL,
  code.medians = NULL,
  marker.pair = NULL,
  use.cluster.counts = T
)
}
\arguments{
\item{dat}{\code{data.table} containing a factored \code{FlowSOM} cluster column}

\item{cluster.col}{Character string; a named column that identifies factored clusters.}

\item{code.medians}{pre-calculated cluster medians; a listed result of \code{fsom.codes.to.clusters}}

\item{marker.pair}{character vector of length 2; if defined, will pre-select the initially displayed dimensions for the bivariate plots.  The default \code{NULL} will pre-select the first two columns of dat for display.}

\item{use.cluster.counts}{Logical. By default, \code{TRUE}; will key \code{dat} by the unique identifier 'sample.id' and cluster column to generate per-sample cluster counts; cluster counts are used as input for factored boxplots.}
}
\value{
\code{shinyapp}
}
\description{
This interactive \code{shiny};\code{shinydashboard} allows a visual exploration of FlowSOM clustering results. \code{plotly} clickable cluster heatmap dimensions and factored datapoints are the primary features.
}
