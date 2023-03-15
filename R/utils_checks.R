sample.id.check<-function(dt.names){
  if(!is.character(dt.names)){
    stop("Need a character vector: 'names(dat)'")
  }
  if(!'sample.id' %in% dt.names){
    stop(paste("Required sample column name not found: 'sample.id' ;",
               "check 'names(dat)' and consider renaming the identifier column with 'data.table::setnames()'",sep="\n"))
  }
}
##
cluster.col.check<-function(dt.names){
  if(!is.character(dt.names)){
    stop("Need a character vector: 'names(dat)'")
  }
  if(!'cluster' %in% dt.names){
    if(any(grepl('cluster',dt.names,ignore.case = T))){
      if(length(which(grepl('cluster',dt.names,ignore.case = T)))==1){
        cluster.col<-grep('cluster',dt.names,ignore.case = T,value=T)
      }else{
        c.names<-paste0(grep('cluster',dt.names,ignore.case = T,value=T),collapse = " ; ")
        stop(paste("More than one cluster column found:",c.names))
      }
    }
  }
  return(cluster.col)
}
##
node.col.check<-function(dt.names){
  if(!'node' %in% dt.names){
    if(any(grepl('node',dt.names,ignore.case = T))){
      if(length(which(grepl('node',dt.names,ignore.case = T)))==1){
        node.col<-grep('node',dt.names,ignore.case = T,value=T)
      }else{
        n.names<-paste0(grep('node',dt.names,ignore.case = T,value=T),collapse = " ; ")
        stop(paste("More than one node column found:",n.names))
      }
    }
  }
  return(node.col)
}
##
