#' @title Cluster FlowSOM codes; calculate per-cluster median expression values
#'
#' @param fsom FlowSOM object, as returned from \code{FlowSOM::SOM}
#' @param k Number of clusters; argument passed to \code{FlowSOM::metaClustering_consensus}
#' @param seed Numerical value; sets random seed; argument passed to \code{FlowSOM::metaClustering_consensus}
#' @param cluster.col Optional character string; if defined, will rename the otherwise generically named 'cluster' column in the median results
#'
#' @details
#' \code{fsom.codes.to.clusters} uses \code{fsom$codes} to first generate metaclusters, then calculates per-cluster median expression values for use primarily in heatmaps/visualizations.
#'
#' @return FlowSOM object appended with list-structured results
#' @export
#'
#'
fsom.codes.to.clusters<-function(fsom,k=20,seed=1337,cluster.col=NULL){
  dims<-colnames(fsom$codes)
  cluster.ls<-list(cluster=list(fac=NULL,dat=list(medians=NULL,medians.scaled=NULL,medians.01=NULL)))
  cluster.ls[['cluster']][['fac']]<-factor(FlowSOM::metaClustering_consensus(fsom$codes, k = 20, seed=seed))
  cluster.ls[['cluster']][['dat']][['medians']]<-data.table::data.table(fsom$codes,cluster=cluster.ls[['cluster']][['fac']])[
    ,lapply(.SD,stats::median),keyby = cluster]
  if(!is.null(cluster.col)) data.table::setnames(cluster.ls[['cluster']][['dat']][['medians']],'cluster',cluster.col)
  cluster.ls[['cluster']][['dat']][['medians.scaled']]<-data.table::copy(cluster.ls[['cluster']][['dat']][['medians']])[
    ,(dims) := lapply(.SD, function(x) (x - mean(x))/stats::sd(x)),.SDcols = dims]
  cluster.ls[['cluster']][['dat']][['medians.01']]<-data.table::copy(cluster.ls[['cluster']][['dat']][['medians']])[
    ,(dims) := lapply(.SD,function(x){scales::rescale(x, from=c(0, max(x)))}),.SDcols = dims]
  return(append(fsom,cluster.ls))
}
