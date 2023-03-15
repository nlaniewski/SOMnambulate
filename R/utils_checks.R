sample.id.check<-function(dt.names){
  if(!'sample.id' %in% dt.names){
    stop(paste("Required sample column name not found: 'sample.id' ;",
               "check 'names(dat)' and consider renaming the identifier column with 'data.table::setnames()'",sep="\n"))
  }
}
