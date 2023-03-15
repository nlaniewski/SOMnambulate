cluster_counts_long <- function(dat,cluster.col=NULL,sample.col=NULL,
                                sample.split.cols=c('subject','visit','condition','batch','batch.date')){
  ##use checks to confirm column names
  sample.id.check(names(dat))
  cluster.col<-cluster.col.check(names(dat))
  #commented lines fail 'devtools::check()' due to 'no visible binding for global variable...'
  #data.table::setkey(dat, sample, cluster)
  data.table::setkeyv(dat, c('sample.id', cluster.col))
  #dat.N.cluster<-dat[data.table::CJ(sample, cluster, unique = TRUE), .N, by = .EACHI]
  dat.N.cluster<-dat[data.table::CJ(get('sample.id'), get(cluster.col), unique = TRUE), .N, by = .EACHI]
  #dat.N.cluster[, prop := N/sum(N) * 100, by = sample]
  dat.N.cluster[, 'prop' := get('N')/sum(get('N')) * 100, by = 'sample.id']
  #dat.N.cluster[, per1million := N*(1E6/sum(N)), by = sample]
  dat.N.cluster[, 'per1million' := get('N')*(1E6/sum(get('N'))), by = 'sample.id']
  ##
  if(!all(sapply(strsplit(dat.N.cluster[['sample.id']],"_"),length)==length(sample.split.cols))){
    stop("Splitting 'sample name' string by '_' results in a length != to length of 'sample.split.cols'")
  }
  dat.N.cluster[,(sample.split.cols) := data.table::tstrsplit(get('sample.id'),"_")]
  return(dat.N.cluster)
}
##
node_counts_long <- function(dat,node.col=NULL,sample.col=NULL,
                             sample.split.cols=c('subject','visit','condition','batch','batch.date')){
  sample.id.check(names(dat))
  node.col<-node.col.check(names(dat))
  #commented lines fail 'devtools::check()' due to 'no visible binding for global variable...'
  #data.table::setkey(dat, sample, node)
  data.table::setkeyv(dat, c('sample.id', node.col))
  #dat.N.node<-dat[data.table::CJ(sample, node, unique = TRUE), .N, by = .EACHI]
  dat.N.node<-dat[data.table::CJ(get('sample.id'), get(node.col), unique = TRUE), .N, by = .EACHI]
  #dat.N.node[, prop := N/sum(N) * 100, by = sample]
  dat.N.node[, 'prop' := get('N')/sum(get('N')) * 100, by = 'sample.id']
  #dat.N.node[, per1million := N*(1E6/sum(N)), by = sample]
  dat.N.node[, 'per1million' := get('N')*(1E6/sum(get('N'))), by = 'sample.id']
  ##
  if(!all(sapply(strsplit(dat.N.node[['sample.id']],"_"),length)==length(sample.split.cols))){
    stop("Splitting 'sample name' string by '_' results in a length != to length of 'sample.split.cols'")
  }
  dat.N.node[,(sample.split.cols) := data.table::tstrsplit(get('sample.id'),"_")]
  return(dat.N.node)
}
