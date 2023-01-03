generate_map_nodes_clusters <- function(dat,dims,scale.dims=T,subsample.val=2E5,k=20){
  codes <- generate_fsom(dat,dims,scale.dims,subsample.val,return.codes = T)
  mc <- generate_clusters(codes,k)
  map <- generate_map_nodes(codes,dat,scale.dims)
  return(data.table::setDT(list(node=map,
                                cluster=mc[map])
                           )
  )
}

generate_fsom <- function(dat,dims,scale.dims=T,subsample.val=2E5,return.codes=F){
  dims<-dims[which(dims %in% names(dat))]
  ##
  if(!is.null(subsample.val)&nrow(dat)>subsample.val){
    set.seed(1337)
    if(scale.dims){
      fsom <- FlowSOM::SOM(as.matrix(dat[sample(.N,subsample.val), lapply(.SD, function(x) (x - mean(x))/stats::sd(x)),.SDcols=dims]),
                           xdim=10,ydim=10
      )
    }else{
      fsom <- FlowSOM::SOM(as.matrix(dat[sample(.N,subsample.val),dims,with=F]),
                           xdim=10,ydim=10
      )
    }
  }else{
    set.seed(1337)
    if(scale.dims){
      fsom <- FlowSOM::SOM(as.matrix(dat[, lapply(.SD, function(x) (x - mean(x))/stats::sd(x)),.SDcols=dims]),
                           xdim=10,ydim=10
      )
    }else{
      fsom <- FlowSOM::SOM(as.matrix(dat[,dims,with=F]),
                           xdim=10,ydim=10
      )
    }
  }
  if(return.codes){
    return(fsom$codes)
  }else{
    return(fsom)
  }
}

generate_clusters <- function(codes,k=20){
  FlowSOM::metaClustering_consensus(codes,k)
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

