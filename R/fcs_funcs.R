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
#' @title Get keyword metadata (!'$P') from .fcs files
#'
#' @param fcs.file.paths Character string; path(s) usually returned from \code{list.files(...,full.names=T,pattern=".fcs")}.
#' @param return.dt Logical. By default, \code{FALSE}; if \code{TRUE}, will return a data.table of keyword/values per .fcs file
#' @param pattern Pattern used to split keyword character strings; specifically for use with barcoded/pooled .fcs files that contain 'collapsed' metadata.
#' @param pattern.unique Logical. By default, \code{TRUE}; returns unique keyword data.tables.
#'
#' @return a list of non-parameter ('$P') keywords per .fcs file; if \code{return.dt=T}, a data.table of keywords/values
#' @export
#'
#'
get.fcs.keywords.metadata <- function(fcs.file.paths,return.dt=F,pattern=NULL,pattern.unique=T){
  fcs.keywords.list <- sapply(flowCore::read.FCSheader(fcs.file.paths),function(h){
    #drop-terms: '$P[0-9]+' (parameters); 'P[0-9]+DISPLAY' (Cytek specific?); Spill(over)
    h<-h[-grep('\\$PAR|\\$P[0-9]+|P[0-9]+DISPLAY|spill',names(h),ignore.case = T)]
    return(as.list(h))
  },simplify = F)
  if(return.dt&is.null(pattern)){
    fcs.keywords.dt <- sapply(fcs.keywords.list, function(kw) data.table::setDT(kw),simplify = F)
  }else if(return.dt&!is.null(pattern)){
    fcs.keywords.dt.split <- lapply(fcs.keywords.list,function(kws){
      lapply(sort(unique(stringr::str_count(kws,pattern))),function(s){
        data.table::as.data.table(sapply(kws[which(stringr::str_count(kws,pattern)==s)],strsplit,pattern))
      })
    })
    if(pattern.unique){
      fcs.keywords.dt.split<-sapply(seq(max(sapply(fcs.keywords.dt.split,length))),function(n){
        fcs.keywords.dt.split<-sapply(fcs.keywords.dt.split,'[',n)
        if(length(unique(fcs.keywords.dt.split))==1) fcs.keywords.dt.split<-fcs.keywords.dt.split[[1]]
        return(fcs.keywords.dt.split)
      })
    }else{
      return(fcs.keywords.dt.split)
    }
  }else{
    return(fcs.keywords.list)
  }
}
#' @title Get a `data.table` of .fcs file paths and related metadata.
#'
#' @param fcs.file.paths Character string; path(s) usually returned from \code{list.files(...,full.names=T,pattern=".fcs")}.
#'
#' @return a `data.table` of full length .fcs file paths and related metadata (parsed from the text header).
#' @export
#'
#'
get.fcs.file.dt<-function(fcs.file.paths){
  f.path<-source.name<-file.size.MB<-NULL
  dt<-data.table::data.table(f.path=fcs.file.paths)
  dt[,source.name:=basename(f.path)]
  dt[,file.size.MB:=signif(file.size(f.path)/1024^2,3)]
  ##
  dts<-get.fcs.keywords.metadata(fcs.file.paths,return.dt = T)
  dt<-cbind(dt,data.table::rbindlist(dts,fill=T))
  #synatically valid names; drop keyword identifier '$'
  names(dt)<-sub("\\$","",names(dt))
  #drop conserved keywords that are non-informative to the user
  drop.terms<-c(
    paste0(c('BEGIN','END'),c('ANALYSIS')),
    paste0(c('BEGIN','END'),c('DATA')),
    paste0(c('BEGINS','ENDS'),c('TEXT')),
    "MODE",
    'BYTEORD',
    'DATATYPE',
    'NEXTDATA',
    'TIMESTEP',
    'APPLY COMPENSATION',
    'USERSETTINGNAME',
    'CHARSET'
  )
  dt<-dt[,!drop.terms[drop.terms %in% names(dt)],with=F]
  #do a few conversions
  for(j in c('DATE')){data.table::set(dt,j=j,value=data.table::as.IDate(dt[[j]],format="%d-%b-%Y"))}
  for(j in c('BTIM','ETIM')){data.table::set(dt,j=j,value=data.table::as.ITime(dt[[j]]))}
  for(j in grep("TOT|DELAY|ASF|VOL$",names(dt),value = T)){data.table::set(dt,j=j,value=as.numeric(dt[[j]]))}
  ##
  return(dt)
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
#' @title reads a .fcs file and converts it to `data.table`
#'
#' @param fcs.dt a `data.table` as returned from `get.fcs.file.dt`
#' @param channel_alias as returned from `get.fcs.channel.alias`
#'
#' @return a `data.table` of raw, un-transformed numeric expression values with character/factor identifier columns
#' @export
#'
#'
fcs.to.dt<-function(fcs.dt,channel_alias=NULL){
  #read .fcs file ('.path')
  #if defined, rename columns using 'channel_alias' as returned from 'get.fcs.channel.alias'
  fcs.tmp<-flowCore::read.FCS(fcs.dt[['f.path']],transformation = F,truncate_max_range = F,
                              channel_alias = if(!is.null(channel_alias)) channel_alias)
  #convert the expression matrix (raw, un-transformed data values) into a data.table
  dt<-data.table::setDT(as.data.frame(fcs.tmp@exprs))
  if(length((fcs.dt[,!'f.path']))>0){
    dt<-cbind(dt,fcs.dt[,!'f.path'])
  }
  return(dt)
}
#' @title a parallelized version of `fcs.to.dt`; essentially a wrapper around `parallel` package functions
#'
#' @param fcs.dt a `data.table` as returned from `get.fcs.file.dt`
#' @param channel_alias as returned from `get.fcs.channel.alias`
#'
#' @return a row-bound `data.table` of raw, un-transformed numeric expression values with character/factor identifier columns
#' @export
#'
#'
fcs.to.dt.parallel<-function(fcs.dt,channel_alias=NULL){
  n.cores<-parallel::detectCores()
  n.paths<-fcs.dt[,.N]
  n<-ifelse(n.paths>n.cores,n.cores,n.paths)
  cl<-parallel::makeCluster(n)
  if(!is.null(channel_alias)){
    parallel::clusterExport(cl,'channel_alias',envir = environment())
  }
  ##
  dt<-data.table::rbindlist(parallel::parLapply(cl,split(fcs.dt,by='f.path'),fcs.to.dt,channel_alias = if(!is.null(channel_alias)) channel_alias))
  ##
  on.exit({parallel::stopCluster(cl);invisible(gc())})
  return(dt)
}
