#' Read a .fcs file and convert to a data.table
#'
#' @param fcs.file.path Character string; 'list.files(...,full.names=T,pattern=".fcs")'
#' @param use.alias Logical. By default, TRUE; will use header-encoded 'marker' names instead of 'channel' names.
#' @param use.alias.split Logical. By default, TRUE; will use 'split' alias names.
#' @param asinh.transform Logical. By default, FALSE; if TRUE, will transform flow/spectral/mass cytometry data using the function 'asinh()'
#' @param cofactor.default Adjust the 'asinh.transform' using a supplied cofactor.  Depending on cytometry data type, use the following default cofactor: flow = 1000, mass = 5
#' @param cofactor.mod Channel/marker-specific adjustment to 'cofactor.default'; provide as a named numeric vector, e.g. 'c("CD3"=2000,"CD8A"=2000)'
#' @param drop.events Logical. By default, TRUE; will drop 'extreme' events using the internal function 'drop_events_count_based()'; these events are most likely to be extreme-negative values introduced due to spectral compensation
#' @param natural.order Logical. By default, TRUE; will sort column names alphabetically then numerically, e.g. "CD3, CD4, CD8, CD45A" instead of "CD3, CD4, CD45A, CD8".
#' @param comp.mat.modified An optional compensation matrix; if provided, will be used in place of a header-encoded matrix ('$SPILL').
#' @param channels.df An optional data.frame returned from the function 'generate.channels.frame()'
#' @param add.identifier.col Logical. By default, TRUE; adds the .fcs filename (basename) as an identifier.
#'
#' @return A data.table composed of .fcs data (cells/events as rows; channels/markers as columns); appended with an additional identifier column ('sample.id').
#' @export
#'
#'
readFCS_dt<-function(fcs.file.path,use.alias=T,use.alias.split=T,
                     asinh.transform=F,cofactor.default=1000,cofactor.mod=NULL,
                     drop.events=T,natural.order=T,comp.mat.modified=NULL,
                     channels.df=NULL,add.identifier.col=T){
  ##
  if(is.null(channels.df)){
    channels.df<-generate.channels.frame(fcs.file.path)
  }else if(!all(c('channels','alias','alias.split') %in% colnames(channels.df))){
    stop(paste("Supplied 'channels.df' should have the following column names: 'channels','alias','alias.split'",
               "Use the returned data.frame from 'generate.channels.frame()'",
               sep = "\n")
    )
  }
  ##
  if(use.alias==F){
    fcs.tmp <- flowCore::read.FCS(fcs.file.path,transformation = F,truncate_max_range = F)
  }else{
    fcs.tmp <- flowCore::read.FCS(fcs.file.path,transformation = F,truncate_max_range = F,
                                  column.pattern = paste0(channels.df$channels,collapse = "|")
    )
  }
  cyto.type<-ifelse(any(grepl("laser",names(fcs.tmp@description),ignore.case = T)),"flow","mass")
  if(cyto.type=="flow"){
    fcs.tmp <- trim.scatter(fcs.tmp)
    if(is.null(comp.mat.modified)){
      fcs.tmp <- flowCore::compensate(fcs.tmp, fcs.tmp@description[[grep("SPILL",names(fcs.tmp@description),ignore.case = T)]])
    }else{
      fcs.tmp <- flowCore::compensate(fcs.tmp,comp.mat.modified)
    }
  }
  ##
  dat <- data.table::as.data.table(fcs.tmp@exprs)
  ##
  if(all(names(dat) == channels.df$channels)){
    if(use.alias==T){
      if(use.alias.split==T){
        data.table::setnames(dat,channels.df$alias.split)
      }else if(use.alias.split==F){
        data.table::setnames(dat,channels.df$alias)
      }
    }
  }
  ##
  if(asinh.transform){
    asinh.cofactors<- stats::setNames(rep(cofactor.default,length(grep("SC|Time|Event_length",names(dat),invert = T))),
                                      nm=grep("SC|Time|Event_length",names(dat),value = T,invert = T)
    )
    if(!is.null(cofactor.mod)){
      if(length(which(names(asinh.cofactors) %in% names(cofactor.mod)))==length(cofactor.mod)){
        for(i in names(cofactor.mod)){
          asinh.cofactors[[i]]<-cofactor.mod[[i]]
        }
      }else{
        stop(paste("Named cofactor(s) not found:"))
      }
    }
    ##
    for (j in names(asinh.cofactors)) data.table::set(dat, j = j, value = asinh(dat[[j]]/asinh.cofactors[[j]]))
  }
  ##
  if(drop.events){
    dat<-dat[!drop_events_count_based(dat)]
  }
  ##
  if(natural.order){
    data.table::setcolorder(dat,stringr::str_sort(names(dat),numeric = T))
  }
  ##
  if(add.identifier.col){
    dat[,'sample.id':=basename(fcs.file.path)]
  }
  return(dat)
}

get.parameters<-function(fcs.file.path){
  header <- flowCore::read.FCSheader(fcs.file.path)[[1]]
  channels <- sapply(c("N","S"),function(i){
    p <- header[grep(paste0("P[0-9]+",i),names(header),value = T)]
    p <- p[order(as.numeric(stringr::str_extract(names(p),"[0-9]+")))]
  },simplify = F)
  ##
  if(any(grepl("[A-Z]{1}[a-z]{1}[0-9]{3}",channels$S))){
    isotope.metal<-function(n){
      metal<-stringr::str_extract(n,"[A-Z]{1}[a-z]{1}")
      isotope<-stringr::str_extract(n,"[0-9]{3}")
      paste0(isotope,metal)
    }
    i<-grep("[A-Z]{1}[a-z]{1}[0-9]{3}",channels$S)
    channels$S[i] <- stringr::str_replace(channels$S[i],"[A-Z]{1}[a-z]{1}[0-9]{3}",isotope.metal)
  }
  ##
  if(any(grepl("140Ce",channels$S))){
    channels$S[grep("140Ce",channels$S)] <- "140Ce_beads"
  }
  ##
  return(channels)
}

generate.channels.frame<-function(fcs.file.path,split.type.position=NULL,dismiss.vars = NULL){
  channels<-get.parameters(fcs.file.path)
  if(is.null(split.type.position)){
    split.type.position<-split.type.position.agnostic(channels$S)
  }
  ##
  channels.retain.vars <- c("Time","Event_length",paste(rep(c("FSC","SSC"),each=3),c("A","H","W"),sep="-"))
  channels.retain <- channels$N[which(channels$N %in% channels.retain.vars)]
  ##
  channels.dismiss.vars <- c("back","noise")#mass cytometers can have named channels that include 'background' and 'noise'
  if(!is.null(dismiss.vars)){
    channels.dismiss.vars <- c(channels.dismiss.vars,dismiss.vars)
  }
  channels.dismiss <- channels$S[grep(paste0(channels.dismiss.vars,collapse = "|"),channels$S,ignore.case = T)]
  ##
  channels.split <- channels$S[which(stringr::str_detect(channels$S,split.type.position$type))]
  channels.split <- channels.split[!channels.split %in% channels.dismiss]
  ##
  channels.keep.n <- sub("S","N",c(names(channels.retain),names(channels.split)))
  ##
  channels.df <- data.frame(channels=channels$N[names(channels$N) %in% channels.keep.n],
                            alias=NA,alias.split=NA)
  channels.df$alias[rownames(channels.df) %in% sub("S","N",names(channels.split))] <- channels.split
  channels.df$alias.split <- sapply(strsplit(channels.df$alias,split.type.position$type),
                                    '[',
                                    split.type.position$position
  )
  ##
  if(any(table(channels.df$alias.split)>1)){
    non.unique<-which(channels.df$alias.split %in% names(which(table(channels.df$alias.split)>1)))
    now.unique<-sapply(channels.df$alias[non.unique],function(i){
      i<-strsplit(i,split.type.position$type)[[1]]
      i<-paste(i[[split.type.position$position]],i[[split.type.position$position-1]],sep='_')
    })
    channels.df$alias.split[non.unique] <- now.unique
  }
  ##
  i<-rownames(channels.df) %in% names(channels.retain)
  channels.df$alias[i] <- channels.df$alias.split[i] <- channels.retain
  ##
  i.to.lower<-grep("intercalator|viability",channels.df$alias.split,ignore.case = T)
  channels.df$alias.split[i.to.lower]<-tolower(channels.df$alias.split[i.to.lower])
  ##
  return(channels.df)
}

split.type.position.agnostic<-function(s){
  common.split.counts <- sapply(c(".", "_", "-", " "),function(split){
    sum(grepl(split,s,fixed = T))
  })
  most.likely.split <- names(common.split.counts)[which.max(common.split.counts)]
  splits<-strsplit(s,most.likely.split)
  split.lengths<-unique(sapply(splits,length))
  if(length(split.lengths)!=1){
    message("Varying split lengths;check naming convention")
  }
  split.position<-which(sapply(seq(split.lengths),function(i){
    any(grepl("CD",sapply(splits,'[',i)))
  }))
  split.type.position<-list(type=most.likely.split,
                            position=split.position
  )
  return(split.type.position)
}

trim.scatter<-function(fcs,trim.low=10000,trim.high=250000){
  scatter.trim.list <- sapply(paste(rep(c("FSC", "SSC"),
                                        each = 3), c("A", "H", "W"), sep = "-"), function(scatter) {
                                          which(fcs@exprs[,scatter]<trim.low|fcs@exprs[,scatter]>trim.high)
                                        })
  fcs@exprs <- fcs@exprs[-Reduce(union, scatter.trim.list),]
  return(fcs)
}

drop_events_median_based<-function(dat,median.cut=0.1){
  if(!data.table::is.data.table(dat)){
    stop("Need data.table (class)")
  }
  transform.these<-grep('SC|Time|Event_length',names(dat),invert = T,value = T,ignore.case = T)
  unique(unlist(dat[, lapply(.SD,function(x){
    b<-seq(min(x),max(x),length.out=301)
    h<-graphics::hist(x,breaks=b,plot=F)
    cut<-stats::median(h$counts[h$counts>0])*median.cut
    cut.index<-which((h$counts>cut))
    list(c(which(x<h$breaks[min(cut.index)]),which(x>h$breaks[max(cut.index)])))
  }), .SDcols=transform.these]))
}

drop_events_count_based<-function(dat){
  if(!data.table::is.data.table(dat)){
    stop("Need data.table (class)")
  }
  transform.these<-grep('SC|Time|Event_length',names(dat),invert = T,value = T,ignore.case = T)
  unique(unlist(dat[, lapply(.SD,function(x){
    b<-seq(min(x),max(x),length.out=301)
    h<-graphics::hist(x,breaks=b,plot=F)
    cut.index<-which((h$counts<=3))
    if(length(cut.index)!=0){
      cut.max.bin<-min(which(sign(h$breaks[cut.index])==1))
      cut.min.bin<-cut.max.bin-1
      cut.min<-h$breaks[cut.index][cut.min.bin]
      cut.max<-h$breaks[cut.index][cut.max.bin]
      list(c(which(x<cut.min),which(x>cut.max)))
    }
  }), .SDcols=transform.these]))
}

channels.frame.check <- function(fcs.file.paths,...){
  channel.frame.lengths<-sapply(fcs.file.paths,function(i){
    suppressMessages(nrow(generate.channels.frame(i)))
  })
  if(length(unique(channel.frame.lengths))>1){
    channels.df <- generate.channels.frame(names(which.max(channel.frame.lengths)),...)
  }else if(all(length(unique(channel.frame.lengths))==1)){
    channels.df <- generate.channels.frame(fcs.file.paths[1],...)
  }
  return(channels.df)
}
