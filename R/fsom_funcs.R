#' @title Cluster FlowSOM codes; calculate per-cluster median expression values
#'
#' @param fsom FlowSOM object, as returned from \code{FlowSOM::SOM}
#' @param k Number of clusters; argument passed to \code{FlowSOM::metaClustering_consensus}
#' @param seed Numerical value; sets random seed; argument passed to \code{FlowSOM::metaClustering_consensus}
#' @param cluster.col Optional character string; if defined, will rename the otherwise generically named 'cluster' column in the median results
#'
#' @details
#' \code{fsom.codes.to.clusters} uses \code{fsom$codes} to first generate metaclusters, then calculates per-cluster median expression values for use primarily in heatmaps/visualizations. If factored clusters already exist (\code{fsom$cluster$fac}), medians will be updated -- useful following \code{fsom.merge.codes}.
#'
#' @return FlowSOM object appended with list-structured results
#' @export
#'
#'
fsom.codes.to.clusters<-function(fsom,k=20,seed=1337,cluster.col=NULL){
  ##
  existing.clusters<-!is.null(fsom$cluster$fac)
  ##
  dims<-colnames(fsom$codes)
  cluster.ls<-list(cluster=list(fac=if(existing.clusters){fsom$cluster$fac}else{NULL},
                                dat=list(medians=NULL,
                                         medians.scaled=NULL,
                                         medians.01=NULL
                                ))
  )
  ##
  if(existing.clusters){
    message("Existing clusters; regenerating medians...")
    cluster.ls$cluster$fac<-fsom$cluster$fac
  }else{
    message(paste("Generating",k,"clusters using FlowSOM::metaClustering_consensus"))
    cluster.ls$cluster$fac <- factor(FlowSOM::metaClustering_consensus(fsom$codes,
                                                                       k = k, seed = seed))
  }
  ##
  cluster.ls$cluster$dat$medians<-data.table::data.table(fsom$codes,cluster=cluster.ls$cluster$fac)[
    ,lapply(.SD,stats::median),keyby = cluster]
  if(!is.null(cluster.col)) data.table::setnames(cluster.ls$cluster$dat$medians,'cluster',cluster.col)
  cluster.ls$cluster$dat$medians.scaled<-data.table::copy(cluster.ls$cluster$dat$medians)[
    ,(dims) := lapply(.SD, function(x) (x - mean(x))/stats::sd(x)),.SDcols = dims]
  cluster.ls$cluster$dat$medians.01<-data.table::copy(cluster.ls$cluster$dat$medians)[
    ,(dims) := lapply(.SD,function(x){scales::rescale(x, from=c(0, max(x)))}),.SDcols = dims]
  ##
  fsom$cluster<-NULL
  return(append(fsom,cluster.ls))
}
##
#' @title Merge FlowSOM codes into a new cluster
#'
#' @param fsom FlowSOM object, as returned from \code{FlowSOM::SOM} and previously clustered using \code{fsom.codes.to.clusters}.
#' @param codes.to.merge Numeric vector; individual codes to be manually merged into a new cluster.
#'
#' @return FlowSOM object with an updated cluster factor (\code{fsom$cluster$fac})
#' @export
#'
#'
fsom.merge.codes<-function (fsom, codes.to.merge)
{
  if (is.null(fsom$cluster$fac)) {
    stop("Need factored clusters; returned from 'fsom.codes.to.clusters(...)")
  }
  m <- fsom$cluster$fac
  levels(m) <- c(levels(m), length(levels(m)) + 1)
  m[codes.to.merge] <- length(levels(m))
  m <- factor(m)
  levels(m) <- c(1:length(levels(m)))
  fsom$cluster$fac<-m
  return(fsom)
}
