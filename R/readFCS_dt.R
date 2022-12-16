readFCS_dt<-function(fcs.file.path,use.alias=T,use.alias.split=T,
                     asinh.transform=F,cofactor.default=1000,cofactor.mod=NULL,
                     drop.events=T,sensible.reorder=T){
  ##
  channels.df<-generate.channels.frame(fcs.file.path)
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
    fcs.tmp <- flowCore::compensate(fcs.tmp,fcs.tmp@description$SPILL)
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
  if(sensible.reorder){
    data.table::setcolorder(dat,cols.sensible.reorder(names(dat)))
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
  channels.df<-data.frame(channels=channels$N,alias=NA,alias.split=NA)
  channels.df$alias[sub("N","S",rownames(channels.df)) %in% names(channels$S)] <- channels$S
  ##
  channels.retain <- c("Time","Event_length")
  channels.df$alias[channels.df$channels %in% channels.retain] <- channels.df$channels[channels.df$channels %in% channels.retain]
  channels.dismiss <- c("back","bead")#"background","bead"
  channels.dismiss.index <- grep(paste0(channels.dismiss,collapse = "|"),channels.df$alias,ignore.case = T)
  channels.df$alias[channels.dismiss.index] <- NA
  channels.df<-channels.df[-which(is.na(channels.df$alias)),]
  ##
  split.these <- which(!channels.df$alias %in% channels.retain)
  channels.df$alias.split[split.these]<-sapply(strsplit(channels.df$alias[split.these],split.type.position$type),
                                               '[',
                                               split.type.position$position
  )
  ##
  non.unique<-which(channels.df$alias.split %in% names(which(table(channels.df$alias.split)>1)))
  now.unique<-sapply(channels.df$alias[non.unique],function(i){
    i<-strsplit(i,split.type.position$type)[[1]]
    i<-paste(i[[split.type.position$position]],i[[split.type.position$position-1]],sep='_')
  })
  channels.df$alias.split[non.unique] <- now.unique
  ##
  channels.df$alias.split[is.na(channels.df$alias.split)]<-channels.df$alias[is.na(channels.df$alias.split)]
  ##
  return(channels.df)
}

split.type.position.agnostic<-function(s){
  common.split.counts <- sapply(c(".", "_", "-", " "),function(split){
    sum(grepl(split,s,fixed = T))
  })
  if(max(common.split.counts)!=length(s)){
    most.likely.split <- ""
  }
  else {
    most.likely.split <- names(common.split.counts)[which.max(common.split.counts)]
  }
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

cols.sensible.reorder<-function(cols){
  scatter<-grep("SC",cols)
  cols[-scatter]<-cols[-scatter][order(stringr::str_extract(cols[-scatter],"[A-Za-z]+"))]
  cols[grep("CD",cols)]<-cols[grep("CD",cols)][order(as.numeric(stringr::str_extract(cols[grep("CD",cols)],"[0-9]+")))]
  cols[grep("IL",cols)]<-cols[grep("IL",cols)][order(as.numeric(stringr::str_extract(cols[grep("IL",cols)],"[0-9]+")))]
  if('Time' %in% cols){
    cols<-cols[-which(cols %in% 'Time')]
    cols<-c(cols,'Time')
  }
  return(cols)
}
