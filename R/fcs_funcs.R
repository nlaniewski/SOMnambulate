#' @title Get parameters ('$P') from .fcs file headers
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
#' @title Get keyword metadata (!'$P') from .fcs file header(s)
#'
#' @param fcs.file.paths Character string; path(s) usually returned from \code{list.files(...,full.names=T,pattern=".fcs")}.
#' @param return.dt Logical. By default, \code{FALSE}; if \code{TRUE}, will return a data.table of keyword/values per .fcs file
#' @param pattern Pattern used to split keyword character strings; specifically for use with barcoded/pooled .fcs files that contain 'collapsed' metadata.
#'
#' @return a list of non-parameter ('$P') keywords per .fcs file; if \code{return.dt=T}, a data.table of keywords/values
#' @export
#'
#'
get.fcs.keywords.metadata <- function(fcs.file.paths,return.dt=F,pattern=NULL){
  fcs.keywords.list <- sapply(flowCore::read.FCSheader(fcs.file.paths),function(h){
    kw<-grep(paste0("\\$",c('B','D','E','M','N','P'),collapse = "|"),names(h),value = T,invert = T)
    return(as.list(h[kw]))
  },simplify = F)
  if(return.dt&is.null(pattern)){
    fcs.keywords.dt <- sapply(fcs.keywords.list, function(kw) data.table::setDT(kw),simplify = F)
  }else if(return.dt&!is.null(pattern)){
    fcs.keywords.dt.split <- lapply(fcs.keywords.list,function(kws){
      lapply(unique(stringr::str_count(kws,pattern)),function(s){
        data.table::as.data.table(sapply(kws[which(stringr::str_count(kws,pattern)==s)],strsplit,pattern))
      })
    })
  }else{
    return(fcs.keywords.list)
  }
}
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
##
#' @title Convert a \code{data.table} of mass cytometry .fcs data into a new .fcs file
#' @description
#' After reading in (usually large amounts of) .fcs data and converting to a \code{data.table}, the process can be reversed and a new .fcs file created; new keywords are written to the header section and parameters updated using \code{FlowCore} functions.
#'
#'
#' @param dt.data A \code{data.table} of .fcs data, as returned from \code{data.table::as.data.table(flowCore::read.FCS(...))}.
#' @param reverse.asinh.cofactor If defined, will reverse the \code{asinh} transformation - using the same cofactor value - typically applied to mass cytometry .fcs data
#' @param keywords.to.add A list of keywords to append to the header section of the .fcs file; if \code{NULL}, returns a 'basic' \code{flowCore::flowFrame}.
#' @param fil If defined, will update the '$FIL' keyword; when writing the new .fcs file, '$FIL' is used to define the \code{flowCore::write.FCS(filename)} argument.
#' @param out.path A \code{file.path}; will write the new .fcs to the defined directory; if \code{NULL}, returns an 'updated' \code{flowCore::flowFrame}.
#'
#' @return Depending on defined arguments, returns a \code{flowCore::flowFrame} with/without updated keywords or writes the .FCS to a defined directory.
#' @export
#'
fcs.from.dt.masscyto<-function(dt.data,reverse.asinh.cofactor=NULL,keywords.to.add=NULL,fil=NULL,out.path=NULL){
  maxRange<-NULL
  if(!is.null(reverse.asinh.cofactor)){
    if(!is.numeric(reverse.asinh.cofactor)) stop("Non-numeric value defined for 'reverse.asinh.cofactor'")
    cols.metal <- grep("[0-9]{3}[A-Z]{1}[a-z]{1}",names(dt.data),value = T)
    dt.data[ , (cols.metal) := lapply(.SD,function(x){sinh(x)*reverse.asinh.cofactor}), .SDcols = cols.metal]
  }
  parms.adf<-data.table::data.table(name = names(dt.data),desc = NA,dt.data[,.(minRange=0,maxRange=unlist(lapply(.SD,function(x) ceiling(max(x)))))])
  parms.adf[,range:=(maxRange+1)]
  ##
  parms.list<-as.list(
    stats::setNames(c(rep(c("32","0,0"),each=nrow(parms.adf)),parms.adf$range),
                    nm=do.call(paste0, expand.grid(paste0('$P',seq(nrow(parms.adf))), c("B","E","R"))))
  )
  ##
  if(!is.null(keywords.to.add)){
    fcs.from.dt<-methods::new("flowFrame",exprs=as.matrix(dt.data),parameters=Biobase::AnnotatedDataFrame(parms.adf),description = parms.list)
    flowCore::keyword(fcs.from.dt)<-c(flowCore::keyword(fcs.from.dt),keywords.to.add)
    flowCore::keyword(fcs.from.dt)[['modified.by']]<-"R_fcs.from.dt.masscyto"
    if(!is.null(fil)){
      flowCore::keyword(fcs.from.dt)[['$FIL']] <- fil
    }
    if(is.null(out.path)){
      return(fcs.from.dt)
    }else{
      if(!dir.exists(out.path)) dir.create(out.path,recursive = T)
      flowCore::write.FCS(fcs.from.dt,file.path(out.path,flowCore::keyword(fcs.from.dt)[['$FIL']]))
    }
  }else{
    return(methods::new("flowFrame",exprs=as.matrix(dt.data),parameters=Biobase::AnnotatedDataFrame(parms.adf),description = parms.list))
  }
}
