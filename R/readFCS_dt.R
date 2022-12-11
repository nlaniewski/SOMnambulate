
readFCS_dt<-function(fcs.file.path,use.alias=T,use.alias.split=T,drop.events=T){
  ##
  channels.df<-generate.channels.frame(fcs.file.path)
  ##
  if(use.alias==F){
    fcs.tmp <- flowCore::read.FCS(fcs.file.path,transformation = F,truncate_max_range = F)
  }else{
    fcs.tmp <- flowCore::read.FCS(fcs.file.path,transformation = F,truncate_max_range = F,channel_alias = channels.df[,c('channels','alias')])
  }
  fcs.tmp <- trim.scatter(fcs.tmp)
  fcs.tmp <- flowCore::compensate(fcs.tmp,fcs.tmp@description$SPILL)
  ##
  dat <- data.table::as.data.table(fcs.tmp@exprs)
  if(use.alias==T){
    if(use.alias.split==T){
      data.table::setnames(dat,channels.df$alias.split)
    }else if(use.alias.split==F){
      data.table::setnames(dat,channels.df$alias)
    }
  }
  ##
  if(drop.events){
    dat<-dat[!drop_events_count_based(dat)]
  }
  ##
  return(dat)
}

generate.channels.frame<-function(fcs.file.path,split.type.position=NULL){
  header <- flowCore::read.FCSheader(fcs.file.path)[[1]]
  channels <- sapply(c("N","S"),function(i){
    p <- header[grep(paste0("P[0-9]+",i),names(header),value = T)]
    p <- p[order(as.numeric(stringr::str_extract(names(p),"[0-9]+")))]
  },simplify = F)
  if(is.null(split.type.position)){
    split.type.position<-split.type.position.agnostic(channels$S)
  }
  channels.df<-data.frame(channels=channels$N,alias=NA)
  channels.df$alias[sub("N","S",rownames(channels.df)) %in% names(channels$S)] <- channels$S
  channels.df$alias[is.na(channels.df$alias)] <- channels.df$channels[is.na(channels.df$alias)]
  channels.df$alias.split<-sapply(strsplit(channels.df$alias,split.type.position$type),'[',split.type.position$position)
  ##
  return(channels.df)
}

split.type.position.agnostic<-function(S){
  common.split.counts <- sapply(c(".", "_", "-", " "),function(split){
    sum(grepl(split,S, fixed = T))
  })
  if(max(common.split.counts)!=length(S)){
    most.likely.split <- ""
  }
  else {
    most.likely.split <- names(common.split.counts)[which.max(common.split.counts)]
  }
  splits<-strsplit(S,most.likely.split)
  split.lengths<-unique(sapply(splits,length))
  if(length(split.lengths)!=1){
    stop("Varying split lengths;check naming convention")
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
  fluors<-grep('SC|Time',colnames(dat),invert = T,value = T,ignore.case = T)
  unique(unlist(dat[, lapply(.SD,function(x){
    b<-seq(min(x),max(x),length.out=301)
    h<-graphics::hist(x,breaks=b,plot=F)
    cut<-stats::median(h$counts[h$counts>0])*median.cut
    cut.index<-which((h$counts>cut))
    list(c(which(x<h$breaks[min(cut.index)]),which(x>h$breaks[max(cut.index)])))
  }), .SDcols=fluors]))
}

drop_events_count_based<-function(dat){
  if(!data.table::is.data.table(dat)){
    stop("Need data.table (class)")
  }
  fluors<-grep('SC|Time',colnames(dat),invert = T,value = T,ignore.case = T)
  unique(unlist(dat[, lapply(.SD,function(x){
    b<-seq(min(x),max(x),length.out=301)
    h<-graphics::hist(x,breaks=b,plot=F)
    cut.index<-which((h$counts<=3))
    cut.max.bin<-min(which(sign(h$breaks[cut.index])==1))
    cut.min.bin<-cut.max.bin-1
    cut.min<-h$breaks[cut.index][cut.min.bin]
    cut.max<-h$breaks[cut.index][cut.max.bin]
    list(c(which(x<cut.min),which(x>cut.max)))
  }), .SDcols=fluors]))
}
