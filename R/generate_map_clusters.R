generate_map_nodes_clusters <- function(dat,dims,scale.dims=T,subsample.val=2E5,k=20){
  codes <- generate_fsom(dat,dims,scale.dims,subsample.val,return.codes = T)
  mc <- generate_clusters(codes,k)
  map <- generate_map_nodes(codes,dat,scale.dims)
  return(data.table::setDT(list(node=map,
                                cluster=mc[map])
                           )
  )
}

#' @title Build a self-organizing map (wrapper)
#' @description A wrapper for \code{FlowSOM::SOM()} with a few \code{data.table} sensibilities for scaling/sampling
#'
#' @param dat Data.table returned from \code{readFCS_dt()}
#' @param dims Character vector of dimensions (data.table names/columns) to use when building the SOMs
#' @param scale.dims Logical. By default, \code{TRUE}; will scale input data before building the SOMs
#' @param subsample.val Optional numeric value defining the number of rows to sample from the input data. The default is \code{NULL}
#' @param ... Additional arguments for \code{FlowSOM::SOM(...)}
#'
#' @return As from \code{FlowSOM::SOM(): "A list containing all parameter settings and results"
#' @export
#'
generate_fsom <- function(dat,dims,scale.dims=T,subsample.val=NULL,...){
  dims<-dims[which(dims %in% names(dat))]
  scale.func<-function(x){(x - mean(x))/stats::sd(x)}
  ##
  fsom<-FlowSOM::SOM(
    data=if(is.null(subsample.val)){
      if(scale.dims){as.matrix(dat[, lapply(.SD,scale.func),.SDcols=dims])
      }else{as.matrix(dat[,dims,with=F])}
    }else{
      if(scale.dims){as.matrix(dat[sample(.N,subsample.val), lapply(.SD,scale.func),.SDcols=dims])
      }else{as.matrix(dat[sample(.N,subsample.val),dims,with=F])}
    },
    ...
  )
}

generate_clusters <- function(codes,k=20,seed=1337){
  FlowSOM::metaClustering_consensus(codes,k,seed)
}

generate_map_nodes <- function(codes,dat,scale.dims=T){
  dims.codes <- colnames(codes)
  if(scale.dims)
    FlowSOM:::MapDataToCodes(codes,as.matrix(dat[, lapply(.SD, function(x) (x - mean(x))/stats::sd(x)),.SDcols=dims.codes]))[,1]
}

cluster.id<-function(dat,cell.type.marker,discriminating.markers){
  if(!cell.type.marker %in% names(dat)){
    stop("Check marker names...")
  }
  if(!all(discriminating.markers %in% names(dat))){
    stop("Check marker names...")
  }
  cluster.medians <- dat[, lapply(.SD, stats::median), keyby = cluster,.SDcols=c(cell.type.marker,discriminating.markers)]
  cid<-which(cluster.medians[,get(cell.type.marker)]>mean(cluster.medians[,get(cell.type.marker)]))
  for(cid.test in discriminating.markers){
    if(length(cid)>1){
      i<-which(cluster.medians[,get(cid.test)]<mean(cluster.medians[,get(cid.test)]))
      if(any(cid %in% i)){
        cid<-cid[which(cid %in% i)]
      }
    }else if(length(cid)==1){
      return(cid)
    }else{
      stop("Need a return value of length 1; check conditionals")
    }
  }
  return(cid)
}

