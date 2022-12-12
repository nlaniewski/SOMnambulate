generate_map_clusters <- function(dat,dims,scale.dims=T,subsample.val=2E5){
  fsom <- generate_fsom()
  mc <- generate_clusters()
  map <- generate_map_nodes <- function(codes=fsom$codes,dat,dims,scale.dims=T){
    if(scale.dims)
      FlowSOM:::MapDataToCodes(codes,as.matrix(dat[, lapply(.SD, function(x) (x - mean(x))/sd(x)),.SDcols=dims]))[,1]
  }
  mc[map]
}

generate_fsom <- function(dat,dims,scale.dims=T,subsample.val=2E5){
  dims<-dims[which(dims %in% names(dat))]
  ##
  if(!is.null(subsample.val)&nrow(dat)>subsample.val){
    set.seed(1337)
    if(scale.dims){
      FlowSOM::SOM(as.matrix(dat[sample(.N,subsample.val), lapply(.SD, function(x) (x - mean(x))/sd(x)),.SDcols=dims]),
                   xdim=10,ydim=10
      )
    }else{
      FlowSOM::SOM(as.matrix(dat[sample(.N,subsample.val),..dims]),
                   xdim=10,ydim=10
      )
    }
  }else{
    set.seed(1337)
    if(scale.dims){
      FlowSOM::SOM(as.matrix(dat[, lapply(.SD, function(x) (x - mean(x))/sd(x)),.SDcols=dims]),
                             xdim=10,ydim=10
      )
    }else{
      FlowSOM::SOM(as.matrix(dat[,..dims]),
                             xdim=10,ydim=10
      )
    }
  }
}

generate_clusters <- function(codes=fsom$codes,k=20){
  FlowSOM::metaClustering_consensus(codes,k)
}

generate_map_nodes <- function(codes=fsom$codes,dat,dims,scale.dims=T){
  if(scale.dims)
  FlowSOM:::MapDataToCodes(codes,as.matrix(dat[, lapply(.SD, function(x) (x - mean(x))/sd(x)),.SDcols=dims]))[,1]
}
