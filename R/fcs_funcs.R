#' Get parameters ('$P') from .fcs file headers
#'
#' @param fcs.file.paths Character string; path(s) usually returned from \code{list.files(...,full.names=T,pattern=".fcs")}.
#' @param return.dt Logical. By default, \code{FALSE}; if \code{TRUE}, will return a data.table of parameters per .fcs file
#'
#' @return a list of parameters per .fcs file; if \code{return.dt=T}, a data.table of parameters
#' @export
#'
#'
get.fcs.parameters<-function(fcs.file.paths,return.dt=F){
  fcs.parameters.list<-sapply(flowCore::read.FCSheader(fcs.file.paths),function(h){
    p.max<-as.numeric(h['$PAR'])
    parm.letters<-unique(sub("\\$P[0-9]+","",grep('\\$P[0-9]',names(h),value = T)))
    p.list<-sapply(parm.letters,function(p){
      pars<-h[grep(paste0("\\$P[0-9]+",p),names(h))]
    },simplify = F)
    for(i in names(which(sapply(p.list,length)<p.max))){
      fill<-setdiff(seq(p.max),as.numeric(stringr::str_extract(names(p.list[[i]]),"[0-9]+")))
      fill.names<-stats::setNames(nm=paste0("$P",fill,i),rep(NA,length(fill)))
      p.list[[i]]<-c(p.list[[i]],fill.names)
    }
    p.list<-lapply(p.list,function(p) p<-p[stringr::str_order(names(p),numeric = T)])
    return(p.list)
  },simplify = F)
  if(return.dt){
    fcs.parameters.list<-sapply(fcs.parameters.list,function(p) data.table::setDT(p),simplify = F)
  }else{
    return(fcs.parameters.list)
  }
}
##
#' @title Get a \code{channel_alias} data.frame from .fcs file headers
#' @description
#' The resultant data.frame is to be used with the \code{channel_alias} argument of \code{flowCore::read.fcs(...)}
#'
#' @param fcs.file.paths Character string; path(s) usually returned from \code{list.files(...,full.names=T,pattern=".fcs")}.
#'
#' @return returns a data.frame containing a 'channels' and 'alias' column
#' @export
#'
#'
get.fcs.channel.alias<-function(fcs.file.paths){
  p.list<-get.fcs.parameters(fcs.file.paths,return.dt = T)
  channel_alias.list<-sapply(p.list,function(dt){
    dt[is.na(S), S := N]
    if(dt[,data.table::uniqueN(S)]!=dt[,.N]){
      stop("Non-unique alias name")
    }else{
      data.table::setnames(dt,c('N','S'),c('channels','alias'))
      return(dt[,.(channels,alias)])
    }
  },simplify = F)
  ##
  if(length(unique(channel_alias.list))==1){
    return(unique(channel_alias.list)[[1]])
  }else{
    warning(paste(
      'channel and/or alias conflict; resolve name discrepancy:',
      paste0(setdiff(
        Reduce(union,lapply(channel_alias.list,'[[','alias')),
        Reduce(intersect,lapply(channel_alias.list,'[[','alias'))),collapse = " : ")
    ),call. = F)
    return(unique(channel_alias.list))
  }
}
