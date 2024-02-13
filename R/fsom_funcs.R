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
#' @param dt A `data.table` as returned from `fcs.to.dt`; the `dt` should be subset to only include numeric columns-of-interest for training the SOM; coerced to matrix.
#' @param .scale Logical: default `FALSE`; if `TRUE`, `dt` will be scaled using an internally-defined function.
#' @param scale.func Default `NULL`; if defined and `.scale=TRUE`, will use the supplied function to override the default internally-defined function.
#' @param seed.val An elite random seed value used to control otherwise random starts.
#' @param ... further arguments passed to `FlowSOM::SOM(...)`
#'
#' @return A list containing `FlowSOM::SOM`-specific parameters/results. Of primary importance is the list element `[['codes']]`.
#' @export
#'
#'
som<-function(dt,.scale=F,scale.func=NULL,seed.val=1337,...){
  if(is.null(scale.func)){
    scale.func<-function(x){(x - mean(x))/stats::sd(x)}
  }
  time.start<-Sys.time()
  set.seed(seed.val)
  fsom<-FlowSOM::SOM(as.matrix(if(.scale){dt[,lapply(.SD,scale.func)]}else{dt}),...,silent=T)
  time.end<-Sys.time()
  elapsed<-as.numeric((time.end)-(time.start))
  if(elapsed<60){tm<-'Seconds'}else{tm<-'Minutes';elapsed<-elapsed/60}
  message(paste(tm,"elapsed:",round(elapsed,3)))
  if(.scale){
    fsom$scale<-TRUE
    fsom$scale.func<-scale.func
  }else{
    fsom$scale<-FALSE
  }
  return(fsom)
}
#' @title Generate a `data.table` of clustered `FlowSOM` codes/SOMs.
#'
#' @param fsom A `FlowSOM::SOM` result as returned from \link{som}.
#' @param append.name Character string; will be used to append the otherwise generically named 'node' and 'cluster' columns.
#' @param k Numeric; if defined, will generate `k` clusters using \link[FlowSOM]{metaClustering_consensus}.
#' @param umap.codes Logical: default `FALSE`; if `TRUE`, the codes/SOMs will be embedded using \link[uwot]{umap}.
#' @param seed.val An elite random seed value used to control otherwise random starts.
#'
#' @return A `data.table` of clustered (factor) codes that represents the direct result of \link[FlowSOM]{SOM}.
#' @export
#'
#'
som.codes.dt<-function(fsom,append.name=NULL,k=NULL,umap.codes=F,seed.val=1337){
  #for R CMD check; data.table vars
  node<-NULL
  #
  if(!"codes" %in% names(fsom)){
    stop("Need a 'FlowSOM' object as returned from 'FlowSOM::BuildSOM()'or 'SOMnambulate::som()'")
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
  return(dt.codes)
}
#' @title Map individual data points to nearest node.
#' @description
#' This function is entirely dependent on the compiled C code found in (hidden function) `FlowSOM:::MapDataToCodes`; as an incredibly useful function, it is used here to map (often tens of millions) data points to their nearest node based on the trained SOM.  The result of `FlowSOM:::MapDataToCodes` is a two column matrix but in this function it is added to `dt` by 'reference' so that no copy is made.
#'
#'
#' @param fsom An object as returned from \link{som}
#' @param dt An object as returned from \link{fcs.to.dt}
#'
#' @return This function modifies `dt` by reference using `data.table`'s `:=` and `set` functions; two derived columns of 'mapped' data -- named generically as 'node' and 'node_dist' -- are added to `dt`. If `fsom` has been list-appended with the return of \link{som.codes.dt}, then an additional column -- 'cluster' -- will be added by reference.
#' @note Assignment operators \link[base]{<-} should not be used with this function as it modifies `dt` by reference; the intent is to avoid copying an often large `dt` object.
#' @references Van Gassen, Sofie, Britt Callebaut, Mary J Van Helden, Bart N Lambrecht, Piet Demeester, Tom Dhaene, and Yvan Saeys. 2015. “FlowSOM: Using Self-Organizing Maps for Visualization and Interpretation of Cytometry Data.” Cytometry Part A 87 (7): 636–45. https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625.
#' @export
#'
#'
map.som.data<-function(fsom,dt){
  #retrieve hidden function from FlowSOM package
  map.data<-utils::getFromNamespace("MapDataToCodes","FlowSOM")
  #retrieve node/cluster names from fsom$codes.dt
  if(!is.null(fsom$codes.dt)){
    node.col<-grep("node",names(fsom$codes.dt),value = T)
    node.cols<-c(node.col,sub("node","node_dist",node.col))
    cluster.col<-grep("cluster",names(fsom$codes.dt),value = T)
  }else{
    node.col<-"node"
    node.cols<-c("node","node_dist")
  }
  #
  message(paste(
    "Mapping data using FlowSOM:::MapDataToCodes;",
    paste("adding columns",paste0("'",node.cols,"'",collapse = " and "),"-- by reference -- to 'dt'."),
    sep = "\n"
  ))
  dt[,(node.cols) := as.list(
    as.data.frame(
      map.data(
        fsom$codes,
        as.matrix(if(fsom$scale){
          dt[,lapply(.SD,fsom$scale.func),.SDcols = colnames(fsom$codes)]
        }else{
          dt[,colnames(fsom$codes),with=F]
        })
      )
    )
  )]
  #
  if(!is.null(fsom$codes.dt)){
    message(paste("Adding column",paste0("'",cluster.col,"'"),"(factor) -- by reference -- to 'dt"))
    data.table::set(dt,j=cluster.col,value = fsom$codes.dt[[cluster.col]][dt[[node.col]]])
  }
  #
  invisible(dt[])
}
