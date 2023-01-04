cluster_counts_long <- function(dat,sample.split.cols=c('subject','visit','condition','batch')){
  if(!all(c('sample','cluster') %in% names(dat))){
    stop("Need 'sample' and 'cluster' columns...")
  }
  data.table::setkey(dat, sample, cluster)
  dat.N.cluster<-dat[data.table::CJ(sample, cluster, unique = TRUE), .N, by = .EACHI]
  dat.N.cluster[, prop := N/sum(N) * 100, by = sample]
  dat.N.cluster[, per1million := N*(1E6/sum(N)), by = sample]
  ##
  if(!all(sapply(strsplit(dat.N.cluster$sample,"_"),length)==length(sample.split.cols))){
    stop("Splitting 'sample' string by '_' results in a length != to length of 'sample.split.cols'")
  }
  dat.N.cluster[,(sample.split.cols) := data.table::tstrsplit(sample,"_")]
  return(dat.N.cluster)
}
##
node_counts_long <- function(dat){
  if(!all(c('sample','node') %in% names(dat))){
    stop("Need 'sample' and 'node' columns...")
  }
  data.table::setkey(dat, sample, node)
  dat.N.node<-dat[data.table::CJ(sample, node, unique = TRUE), .N, by = .EACHI]
  dat.N.node[, prop := N/sum(N) * 100, by = sample]
  dat.N.node[, per1million := N*(1E6/sum(N)), by = sample]
  ##
  if(!all(sapply(strsplit(dat.N.node$sample,"_"),length)==length(sample.split.cols))){
    stop("Splitting 'sample' string by '_' results in a length != to length of 'sample.split.cols'")
  }
  dat.N.node[,(sample.split.cols) := data.table::tstrsplit(sample,"_")]
  return(dat.N.node)
}
