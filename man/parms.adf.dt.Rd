% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_funcs.R
\name{parms.adf.dt}
\alias{parms.adf.dt}
\title{Generate a \code{data.frame} of .fcs parameters from a \code{data.table} of .fcs expression data}
\usage{
parms.adf.dt(dt, to.numeric = NULL, name.fix = NULL)
}
\arguments{
\item{dt}{a \code{data.table} as returned from \link{fcs.to.dt}.}

\item{to.numeric}{Character vector; will attempt to convert each string element (named column in \code{dt}) to numeric.}

\item{name.fix}{a \code{data.table} with two columns: 'name' and 'name.fix'; if 'name' matches, the value will be replaced with 'name.fix'.}
}
\value{
a \code{data.frame} of .fcs parameters; used to define a (new) \link[flowCore]{flowFrame}'s parameter slot using \link[Biobase]{AnnotatedDataFrame}.
}
\description{
As a utility function, this returns a \code{data.frame} of .fcs parameters for use in building a new .fcs file by fulfilling the \link[Biobase]{AnnotatedDataFrame} slot. Only numeric columns will be included.
}
\examples{
#from the 'get.fcs.file.dt' example:
#from the 'get.fcs.channel.alias' example:
#from the 'fcs.to.dt' example:

dt<-SOMnambulate:::prepared.examples('dt')

parms.adf<-SOMnambulate:::parms.adf.dt(dt,to.numeric='aliquot.seq')
parms.adf[]

}
