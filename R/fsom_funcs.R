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
#' @title Build a self-organized map (SOM)
#'
#' @description
#' Essentially a wrapper for \link[FlowSOM:SOM]{FlowSOM::SOM} with a few \link[data.table]{data.table}-specific operations.
#'
#' @param dt A \link[data.table]{data.table} as returned from \link{fcs.to.dt}; the `dt` should be subset to only include numeric columns-of-interest for training the SOM; coerced to matrix.
#' @param .scale Logical: default `FALSE`; if `TRUE`, `dt` will be scaled using an internally-defined function.
#' @param seed.val An elite random seed value used to control otherwise random starts.
#' @param ... Further arguments passed to \link[FlowSOM:SOM]{FlowSOM::SOM}
#'
#' @return A list containing \link[FlowSOM:SOM]{FlowSOM::SOM}-specific parameters/results. Of primary importance is the list element `[['codes']]`.
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples(example.type='dt')
#'
#' pbmc.markers<-c("CD3","CD4","CD8a","CD14","CD19","CD56","TCRgd")
#' pbmc.dims<-grep(paste0(pbmc.markers,"_",collapse="|"),names(dt),value = TRUE)
#' for(j in pbmc.dims){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#'
#' fsom<-som(dt[,pbmc.dims,with=FALSE],.scale=TRUE,map=FALSE)
#' str(fsom)
#'
som<-function(dt,.scale=FALSE,seed.val=1337,...){
  if(.scale){
    scale.func<-function(x){(x - mean(x))/stats::sd(x)}
  }
  ts<-Sys.time()
  set.seed(seed.val)
  message("Building a self-organizing map using FlowSOM::SOM...")
  fsom<-FlowSOM::SOM(as.matrix(if(.scale){dt[,lapply(.SD,scale.func)]}else{dt}),...,silent=T)
  te<-Sys.time()
  message(paste("Minutes elapsed:",round(difftime(te, ts, units = "mins"),2)))
  if(.scale){
    fsom$scale<-TRUE
  }else{
    fsom$scale<-FALSE
  }
  return(fsom)
}
#' @title Generate a \link[data.table]{data.table} of clustered \link[FlowSOM:SOM]{FlowSOM::SOM} codes/SOMs.
#' @description
#' Using the direct result of \link{som}, the `[['codes']]` are clustered and a column-bound `data.table` is generated.  As it is computationally expensive to UMAP-embed often millions of data points, an alternative is to embed the codes themselves for added interpretation of clustering results.
#'
#' @param fsom A \link[FlowSOM:SOM]{FlowSOM::SOM} result as returned from \link{som}.
#' @param append.name Character string; will be used to append the otherwise generically named 'node' and 'cluster' columns.
#' @param k Numeric; if defined, will generate `k` clusters using \link[ConsensusClusterPlus]{ConsensusClusterPlus}.
#' @param umap.codes Logical: default `FALSE`; if `TRUE`, the codes/SOMs will be embedded using \link[uwot]{umap}.
#' @param seed.val An elite random seed value used to control otherwise random starts.
#' @param pItem Argument passed to \link[ConsensusClusterPlus]{ConsensusClusterPlus}; directly influences final clustering result.
#' @param reps Argument passed to \link[ConsensusClusterPlus]{ConsensusClusterPlus}; directly influences final clustering result.
#'
#' @return An updated `fsom`; the new list element `[['codes.dt']]` will contain a \link[data.table]{data.table} of clustered (factor) codes using the direct result of \link[FlowSOM]{SOM}.
#' @export
#'
#' @examples
#' fsom<-SOMnambulate:::prepared.examples(example.type='fsom')
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=10,umap.codes = TRUE)
#'
#' head(fsom$codes.dt[])
#'
#' ggplot2::ggplot(fsom$codes.dt,ggplot2::aes(umap.1,umap.2)) +
#' ggplot2::geom_point(ggplot2::aes(color=cluster_pbmc))
#'
som.codes.dt<-function(fsom,append.name=NULL,k=NULL,umap.codes=FALSE,seed.val=1337,pItem=1,reps=100){
  #
  if(!"codes" %in% names(fsom)){
    stop("Need a 'FlowSOM' object as returned from 'FlowSOM::BuildSOM()'or 'SOMnambulate::som()'")
  }
  dt.codes<-data.table::data.table(fsom$codes)
  dt.codes[,node:=seq(.N)]
  if(!is.null(k)){
    message(paste("Generating", k, "clusters using ConsensusClusterPlus::ConsensusClusterPlus(...)"))
    dt.codes[,cluster:=factor(
      suppressMessages(
        ConsensusClusterPlus::ConsensusClusterPlus(
          d = t(fsom$codes),
          maxK = k,
          pItem = pItem,
          reps = reps,
          title=tempdir(),
          plot = "pdf",
          distance = "euclidean",
          seed=seed.val,
          verbose = FALSE
        ))[[k]]$consensusClass
      )]
    #dt.codes[,cluster:=factor(FlowSOM::metaClustering_consensus(fsom$codes,k=k,seed=seed.val))]
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
  fsom$codes.dt<-dt.codes
  return(fsom)
}
#' @title Map individual data points to nearest node.
#' @description
#' This function is entirely dependent on the compiled C code found in (hidden function) `FlowSOM:::MapDataToCodes`; it is used here to map (often tens of millions) data points to their nearest node based on a trained SOM.  The result of `FlowSOM:::MapDataToCodes` is a two column matrix but in this function it is added (column-bound) to `dt` by 'reference' so that no copy (in memory) is made.
#'
#' @param fsom An object as returned from \link{som} or \link{som.codes.dt}
#' @param dt An object as returned from \link{fcs.to.dt}
#'
#' @return This function modifies `dt` by reference using \link[data.table]{data.table}'s \link[data.table]{:=} and \link[data.table]{set} functions; two derived columns of 'mapped' data -- named generically as 'node' and 'node_dist' -- are added to `dt`. If `fsom` has been list-appended with the return of \link{som.codes.dt}, then an additional column -- 'cluster' -- will be added by reference.
#' @note Assignment operators \link[base]{<-} should not be used with this function as doing so will make a copy of the often large (in-memory) `dt`; instead, `dt` is modified/updated by reference.
#' @references Van Gassen, Sofie, Britt Callebaut, Mary J Van Helden, Bart N Lambrecht, Piet Demeester, Tom Dhaene, and Yvan Saeys. 2015. “FlowSOM: Using Self-Organizing Maps for Visualization and Interpretation of Cytometry Data.” Cytometry Part A 87 (7): 636–45. https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625.
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples(example.type='dt')
#'
#' pbmc.markers<-c("CD3","CD4","CD8a","CD14","CD19","CD56","TCRgd")
#' pbmc.dims<-grep(paste0(pbmc.markers,"_",collapse="|"),names(dt),value = TRUE)
#' for(j in pbmc.dims){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#'
#' fsom<-SOMnambulate:::prepared.examples(example.type='fsom')
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=10,umap.codes = TRUE)
#'
#' #'dt' as returned from 'fcs.to.dt()'; expression values and metadata
#' head(dt[])
#'
#' #do not assign to environment; the function updates 'dt' by reference
#' map.som.data(fsom,dt)
#'
#' #three new columns added -- by reference -- to 'dt'
#' head(dt[])
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
  if(fsom$scale){
    scale.func<-function(x){(x - mean(x))/stats::sd(x)}
  }
  dt[,(node.cols) := as.list(
    as.data.frame(
      map.data(
        fsom$codes,
        as.matrix(if(fsom$scale){
          dt[,lapply(.SD,scale.func),.SDcols = colnames(fsom$codes)]
        }else{
          dt[,colnames(fsom$codes),with=F]
        })
      )
    )
  )]
  #
  if(!is.null(fsom$codes.dt)&length(cluster.col)==1){
    message(paste("Adding column",paste0("'",cluster.col,"'"),"(factor) -- by reference -- to 'dt"))
    data.table::set(dt,j=cluster.col,value = fsom$codes.dt[[cluster.col]][dt[[node.col]]])
  }
  #
  invisible(dt[])
}
#' @title Merge two or more \link{som.codes.dt} clusters
#' @description
#' Due to over-clustering, some resultant clusters may need to be merged. The decision to merge clusters is based on multiple factors (including but not limited to): heatmap expression values of shared features (dendrogram), primary data (scatterplots), and population-level associations (cluster counts).
#'
#' @param fsom An object as returned from \link{som.codes.dt}.
#' @param clusters.to.merge A list; each list element (numeric vector) must be a two or more clusters to merge.
#' @param merge.as.nodes Logical: default `FALSE`; if `TRUE`, the supplied `clusters.to.merge` will merge nodes instead of clusters; consider this a special-use case.
#' @param preserve.factor.levels Logical: default `TRUE`; preserves existing factor levels.
#'
#' @return A modified fsom object; the cluster column in `fsom$codes.dt` will be updated -- by reference -- with the merged results. This function modifies `fsom$codes.dt` by reference using \link[data.table]{data.table}'s \link[data.table]{:=} and \link[data.table]{set} functions.
#' @note Assignment operators \link[base]{<-} should not be used with this function; instead, `fsom$codes.dt` is modified/updated by reference.
#' @export
#'
#' @examples
#' fsom<-SOMnambulate:::prepared.examples(example.type='fsom')
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=10,umap.codes = FALSE)
#' fsom.dup<-data.table::copy(fsom)
#'
#' #current clusters (factor)
#' fsom$codes.dt[['cluster_pbmc']]
#'
#' #for the sake of this example, assume the following clusters need to be merged:
#' clusters.to.merge<-list(c(1,2),c(9,10))
#' merge_som.clusters(fsom,clusters.to.merge)
#' merge_som.clusters(fsom.dup,clusters.to.merge,preserve.factor.levels=FALSE)
#'
#' #result of 'preserve.factor.levels' argument
#' levels(fsom$codes.dt[['cluster_pbmc']])
#' levels(fsom.dup$codes.dt[['cluster_pbmc']])
#'
#' #updated clusters (factor)
#' fsom$codes.dt[['cluster_pbmc']]
#'
#' #special-use case; merging individual nodes
#' merge_som.clusters(fsom,list(c(1,2,3)),merge.as.nodes=TRUE)
#' fsom$codes.dt[['cluster_pbmc']]
#'
merge_som.clusters<-function(fsom,clusters.to.merge,merge.as.nodes=F,preserve.factor.levels=TRUE){
  if(is.null(fsom$codes.dt)){
    stop("Need fsom as returned from som.codes.dt()")
  }else{
    cluster.col<-grep("cluster",names(fsom$codes.dt),value = T)
  }
  if(!is.list(clusters.to.merge)){
    stop(paste(
      "clusters.to.merge needs to be a list;",
      "each list element (numeric vector of cluster numbers) will be merged into a new cluster",
      sep = "\n")
    )
  }
  ##use cluster values in clusters.to.merge to generate an index of which nodes to re-factor
  if(!merge.as.nodes){
    nodes.to.merge<-lapply(clusters.to.merge,function(i){
      fsom$codes.dt[,.I[get(cluster.col) %in% i]]
    })
  }else{
    nodes.to.merge<-clusters.to.merge
  }
  ##use nodes.to.merge indices to re-factor and re-level existing clusters
  for(nodes in nodes.to.merge){
    m <- fsom$codes.dt[[cluster.col]]
    levels(m) <- c(levels(m), length(levels(m)) + 1)
    m[nodes] <- length(levels(m))
    m <- factor(m)
    if(!preserve.factor.levels){
      levels(m) <- c(1:length(levels(m)))
    }
    ##
    data.table::set(fsom$codes.dt,j=cluster.col,value = m)
  }
  ##
  message(paste("Merging clusters and updating",paste0("'",cluster.col,"'"), "-- by reference -- to fsom$codes.dt"))
}
