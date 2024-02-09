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
#' @title Build a self-organized map (SOM)
#'
#' @description
#' Essentially a wrapper for `FlowSOM::SOM` with a few `data.table`-specific operations.
#'
#'
#' @param dt A `data.table` as returned from `fcs.to.dt`; the `dt` should be subset to only include columns-of-interest for training the SOM; coerced to matrix.
#' @param .scale Logical: default `FALSE`; if `TRUE`, the (subset) `dt` will be scaled using an internally-defined function.
#' @param scale.func Default `NULL`; if defined and `.scale=TRUE`, will use the supplied function to override the default internally-defined function.
#' @param seed.val An elite random seed value used to control otherwise random starts.
#' @param ... further arguments passed to `FlowSOM::SOM(...)`
#'
#' @return A `FlowSOM` object; a list containing parameters/results.
#' @export
#'
#'
som<-function(dt,.scale=F,scale.func=NULL,seed.val=1337,...){
  if(is.null(scale.func)){
    scale.func<-function(x){(x - mean(x))/stats::sd(x)}
  }
  time.start<-Sys.time()
  set.seed(seed.val)
  fsom<-FlowSOM::SOM(as.matrix(if(.scale){dt[,lapply(.SD,scale.func)]}else{dt}),...)
  time.end<-Sys.time()
  elapsed<-as.numeric((time.end)-(time.start))
  if(elapsed<60){tm<-'Seconds'}else{tm<-'Minutes';elapsed<-elapsed/60}
  message(paste(tm,"elapsed:",round(elapsed,3)))
  class(fsom) <- c(class(fsom),"FlowSOM")
  return(fsom)
}
#' @title Update a `FlowSOM` object with a `data.table` of clustered codes/SOMs.
#'
#' @param fsom A `FlowSOM` object as returned from `som`.
#' @param append.name Character string; will be used to append the otherwise generically named 'node' and 'cluster' columns.
#' @param k Numeric; if defined, will generate `k` clusters using `FlowSOM::metaClustering_consensus()`.
#' @param umap.codes Logical: default `FALSE`; if `TRUE`, the codes/SOMs will be embedded using `uwot::umap()`.
#' @param seed.val An elite random seed value used to control otherwise random starts.
#'
#' @return A list appended `FlowSOM`object updated with a `$code.dt` containing clusters.
#' @export
#'
#'
fsom.codes.dt<-function(fsom,append.name=NULL,k=NULL,umap.codes=F,seed.val=1337){
  #for R CMD check; data.table vars
  node<-NULL
  #
  if(!"FlowSOM" %in% class(fsom)){
    stop("Need a 'FlowSOM' object as returned from 'FlowSOM::BuildSOM(())'or 'SOMnambulate::som'")
  }
  dt.codes<-data.table::data.table(fsom$codes)
  dt.codes[,node:=seq(.N)]
  if(!is.null(k)){
    message(paste("Generating", k, "clusters using FlowSOM::metaClustering_consensus(...)"))
    dt.codes[,cluster:=factor(FlowSOM::metaClustering_consensus(fsom$codes,k=k,seed=seed.val))]
  }
  if(umap.codes){
    message("UMAP embedding of 'fsom$codes' using uwot::umap(...)")
    set.seed(seed.val)
    dt.codes[,paste0('umap.',1:2) := as.list(as.data.frame(uwot::umap(fsom$codes)))]
  }
  if(!is.null(append.name)){
    old<-grep("node|cluster",names(dt.codes),value = T)
    new<-paste(old,append.name,sep="_")
    data.table::setnames(dt.codes,old,new)
  }
  fsom<-append(fsom,list(codes.dt=dt.codes))
  return(fsom)
}
