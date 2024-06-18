#' @title Get local maxima and minima from a density distribution
#' @description
#' Get local maxima (peaks) and minima (valleys) from a density distribution
#'
#'
#' @param x Numeric vector.
#' @param plot Logical; if `TRUE`, a plot will be drawn displaying a \link[stats]{density} distribution with vertical \link[graphics]{abline}s: "purple peaks" and "violet valleys".
#' @param quantile.probs Numeric vector of length 2; used to define lower/upper quantile \link[stats:quantile]{probabilities}; trims 'extreme events' before calculating the density distribution.
#'
#' @return a list; contains numeric values: local maxima `[['lmaxima']]` and local minima `[['lminima']]`.
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples('dt')
#'
#' dt[,CD45_112Cd := asinh(CD45_112Cd/10)]
#' x<-dt[['CD45_112Cd']]
#'
#' #maxima and minima values returned for extreme tails
#' SOMnambulate:::get.local.maxima.minima(x,plot=TRUE)
#'
#' #quantile 'trim'
#' SOMnambulate:::get.local.maxima.minima(x,quantile.probs=c(0,0.999),plot=TRUE)
get.local.maxima.minima<-function(x,quantile.probs=NULL,plot=F){
  d <- stats::density(
    if(!is.null(quantile.probs)){
      q<-stats::quantile(x,quantile.probs)
      x[x>q[1] & x<q[2]]
    }else{
      x
    }
  )
  diff1<-sign(diff(d$y))
  diff2<-diff(diff1)
  #
  lmaxima<-which(diff2==(-2))+1
  lminima<-which(diff2==(2))+1
  #
  lmaxima.val<-d$x[lmaxima]
  lminima.val<-d$x[lminima]
  #
  if(plot){
    plot(d,main="Density Distribution")
    graphics::abline(v=lmaxima.val,col="purple",lwd=2)
    graphics::abline(v=lminima.val,col="violet",lwd=2)
  }
  return(list(lmaxima=lmaxima.val,
              lminima=lminima.val))
}
#' @title Generate a barcode key
#' @description
#' Based on a 'x-choose-m' scheme, a full length barcode key -- all possible combinations -- is generated. For mass cytometry live-cell barcoding, the scheme is usually '6-choose-3' (20 unique combinations) or '7-choose-3' (35 unique combinations).
#'
#' @param barcode.dims Character vector.
#' @param m Argument passed to \link[utils]{combn}.
#'
#' @return A barcode key (\link[data.table]{data.table}).
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples("dt")
#'
#' barcode.dims<-grep("CD45_",names(dt),value=TRUE)
#'
#' SOMnambulate:::key.barcode(barcode.dims)[]
key.barcode<-function(barcode.dims,m=3){
  .length<-length(barcode.dims)
  key<-t(utils::combn(length(barcode.dims),m))
  .nrow<-nrow(key)
  #
  barcode.key<-matrix(0,nrow=.nrow,ncol=.length,dimnames=list(NULL,barcode.dims))
  #
  for(i in seq(.nrow)){
    barcode.key[i,key[i,]]<-1
  }
  #
  barcode.key<-data.table::data.table(barcode.key)
  barcode.key[,barcode:=seq(.N)]
  #
  invisible(barcode.key)
}
#' @title Generate a 'keyed' `[['codes']]`
#' @description
#' Based on the 'x-choose-m' scheme used to barcode individual samples, this function will generate a unique barcode value for each SOM; the barcode value can then be assigned to the primary input data allowing for identification of sample-specific events/cells.
#'
#' @param fsom Object; the result of \link{som.codes.dt}.  The `[['codes']]` result will be used to generate a 'code-specific' barcode key.
#' @param m Argument passed to \link{key.barcode}.
#' @param quantile.probs Argument passed to \link{get.local.maxima.minima}; Numeric vector of length 2; used to define lower/upper quantile probabilities; trims 'extreme events' before calculating the density distribution.
#'
#' @return A merged \link[data.table]{data.table}; the 'keyed' `[['codes']]` are merged with the result of \link{key.barcode}. Of primary importance is the 'barcode' column which is used to assign barcode values to cellular events (primary data).
#' @export
#'
#' @examples
#'
#' dt<-SOMnambulate:::prepared.examples("dt")
#'
#' barcode.dims<-grep("CD45_",names(dt),value=TRUE)
#' for(j in barcode.dims){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#' fsom<-som(dt[,barcode.dims,with=FALSE],map=FALSE)
#' fsom<-som.codes.dt(fsom,append.name = 'barcode')
#'
#' barcodes<-key.codes(fsom)
#' barcodes[]
#'
#' SOMnambulate::map.som.data(fsom,dt)
#' dt[,barcode := barcodes$barcode[node_barcode]]
#' dt[,.(node_barcode,barcode)]
#' dt[,sort(unique(barcode))]
#' dt[,.N,keyby=.(barcode,stim.condition)]
#'
key.codes<-function(fsom,m=3,quantile.probs=NULL){
  .barcode.dims<-colnames(fsom$codes)
  .nrow<-nrow(fsom$codes)
  .length<-length(.barcode.dims)
  barcode.dt<-data.table::data.table(matrix(0,nrow=.nrow,ncol=.length,dimnames = list(NULL,.barcode.dims)))
  #
  pv<-sapply(.barcode.dims,function(j){
    pv<-get.local.maxima.minima(fsom$codes.dt[[j]],quantile.probs)
    #assumption: near-zero/zero peak will always be the first element
    pv$lmaxima <- pv$lmaxima[1]
    #assumption: valley-of-interest (that which follows the near-zero/zero) will always be the first element
    pv$lminima <- pv$lminima[1]
    #
    pv$half.low<-pv$lminima-(pv$lminima-pv$lmaxima)/2
    return(pv)
  },simplify = F)
  #
  for(j in .barcode.dims){
    data.table::set(
      barcode.dt,
      i=fsom$codes.dt[,.I[get(j)>pv[[j]]$lminima]],
      j=j,
      value=1)
  }
  #
  barcode.dt[,barcode.N := rowSums(.SD),.SDcols = .barcode.dims]
  #reassignment
  for(i in barcode.dt[,.I[barcode.N==2]]){
    for(j in .barcode.dims){
      data.table::set(
        barcode.dt,
        i=i,
        j=j,
        value=ifelse(fsom$codes.dt[i,get(j)]>pv[[j]]$half.low,1,0))
    }
  }
  #
  barcode.dt[,barcode.N := rowSums(.SD),.SDcols = .barcode.dims]
  #
  barcode.dt<-merge(barcode.dt,key.barcode(.barcode.dims,m),all.x=T,sort=F)
  barcode.dt[is.na(barcode),barcode:=0]
  #
  invisible(barcode.dt)
}
#' @title Generate a barcode key
#'
#' @param barcode.dims Character vector
#' @param m argument passed to \code{utils::combn(x,m)}
#'
#' @return a barcode key (matrix)
#'
#'
#'
get.barcode.key<-function(barcode.dims,m=3){
  n <- length(barcode.dims)
  barcode.combinations <- t(utils::combn(n, m));l<-nrow(barcode.combinations)
  barcode.key <- matrix(data = 0, nrow = l,
                        ncol = n, dimnames = list(NULL,barcode.dims))
  for (i in seq(l)) {
    barcode.key[i,barcode.combinations[i,]]<-1
  }
  return(barcode.key)
}
#' @title Find a density distribution 'valley'
#' @description
#' Using the second derivative of the differences in a density distribution's y values, both peaks and valleys can be detected; if a single valley is detected, its x value is returned. A standard deviation adjustment is included to 'adjust/shift' the valley -- used during 'debarcoding' to fine-tune assignment results.
#'
#'
#' @param x Numeric vector
#' @param quantile.trim Logical. By default, \code{TRUE}; trims 'extreme' values.
#' @param no.valley.override Numeric. By default, \code{NULL}; If defined and no valley is detected, this value will be used in place.
#'
#' @return Numeric representing the 'valley' value in a density distribution.
#'
#'
#'
get.valley<-function(x,quantile.trim=T,no.valley.override=NULL){
  if(quantile.trim){
    q.vals <- stats::quantile(x, probs = c(0.001, 0.999))#trim 'extreme' values
    d <- stats::density(x[x > q.vals[1] & x < q.vals[2]])
  }else{
    d <- stats::density(x)
  }
  diff.derivative.1<-sign(diff(d$y))#(sign of) derivative 1 of density 'y' values; monotonic 'slopes'
  diff.derivative.2<-diff(diff.derivative.1)#derivative 2; peaks and valleys
  valley<-which(diff.derivative.2==2)#-2 == peak(s);2 == valley
  if(length(valley)!=1){
    valley<-valley[which.min(d$y[valley])]
  }
  if(length(valley)==1){
    return(d$x[valley])
  }else{
    if(!is.null(no.valley.override)){
      return(no.valley.override)
    }else{
      stop("Can't get a valley value; check density distribution")
    }
  }
  # if(length(valley)==1){
  #   if(!is.null(sd.adjust)){
  #     return(d$x[valley]-stats::sd(d$y)*sd.adjust)#return x value where valley occurs minus 'sd.adjust' multiples of y value standard deviations; 'expression' value
  #   }else{
  #     return(d$x[valley])#return x value where valley occurs; 'expression' value
  #   }
  # }else{
  #   stop(paste(paste(length(valley),"Valleys detected..."),
  #              paste("Valleys at:",paste(round(d$x[valley],5),collapse = ' ; ')),
  #              sep='\n')
  #   )
  # }
}
#' @title Find peaks and valley(s) from a density distribution
#'
#' @param x Numeric vector
#' @param quantile.trim Logical. By default, \code{TRUE}; trims 'extreme' values.
#' @param sd.adjust Numeric; if defined, will adjust the return value by multiples of density value standard deviations.
#'
#' @return a named list: \code{$peaks} (peak values) and \code{$valley} (valley values)
#'
#'
#'
get.peaks.and.valley<-function(x,quantile.trim=T,sd.adjust=NULL){
  if(quantile.trim){
    q.vals <- stats::quantile(x, probs = c(0.001, 0.999))#trim 'extreme' values
    d <- stats::density(x[x > q.vals[1] & x < q.vals[2]])
  }else{
    d <- stats::density(x)
  }
  diff.derivative.1<-sign(diff(d$y))#(sign of) derivative 1 of density 'y' values; monotonic 'slopes'
  diff.derivative.2<-diff(diff.derivative.1)#derivative 2; peaks and valleys
  ##
  peaks<-which(diff.derivative.2==(-2))#derivative 2 peaks are at -2
  peaks.val<-d$x[peaks]
  ##
  valley<-which(diff.derivative.2==2)#derivative 2 valley(s) are at 2
  ##
  if(length(valley)==1){
    if(!is.null(sd.adjust)){
      valley.val <- d$x[valley]-stats::sd(d$y)*sd.adjust#return x value where valley occurs minus 'sd.adjust' multiples of y value standard deviations; 'expression' value
    }else{
      valley.val <- d$x[valley]#return x value where valley occurs; 'expression' value
    }
  }else{
    message(paste(paste(length(valley),"Valleys detected..."),
                  paste("Valleys at:",paste(round(d$x[valley],5),collapse = ' ; ')),
                  sep='\n')
    )
    valley.val<-d$x[valley]
  }
  return(list(peaks=peaks.val,
              valley=valley.val)
  )
}
#' @title Assign barcode values to the 'codes' result of a \code{FlowSOM::SOM} object based on an internally generated 'key'
#' @description
#' Based on the dimensions/markers used to barcode individual samples, this function will generate a unique barcode value for each SOM; the barcode value can then be assigned to the primary input data allowing for identification of sample-specific events/cells.
#'
#'
#' @param fsom Object; the result of \code{FlowSOM::SOM(...)}.  The 'codes' (SOMs) result will be used to generate 'code-specific' barcode values.
#' @param m Numeric. Default = 3. Argument passed to \code{utils::combn(x,m)}; used to define the barcode 'scheme': x-choose-m
#'
#' @return \code{fsom} object -- modified; list element (\code{fsom$barcodes}) with the barcoding result (numeric vector).
#' @export
#'
#'
barcode.assignment.fsom.codes<-function(fsom,m=3){
  n<-ncol(fsom$codes)
  barcode.key<-get.barcode.key(colnames(fsom$codes),m=m)
  valleys <- apply(fsom$codes,2,get.valley)
  node.key <- sapply(names(valleys), function(x) {
    key <- ifelse(fsom$codes[, x] > valleys[x], 1, 0)
  })
  node.key.t <- t(node.key)
  barcode.key.t<-t(barcode.key)
  barcode.assignment <- apply(node.key.t, 2, function(i) {
    if (max(colSums(i == barcode.key.t)) == n) {
      barcode <- which(colSums(i == barcode.key.t) == n)
    }
    else {
      barcode <- 0
    }
  })
  fsom$barcodes<-barcode.assignment
  return(fsom)
}
#' @title Generate barcode assignments based on the result of \code{FlowSOM::SOM}
#'
#' @param dat A \code{data.table} of transformed .fcs data containing the markers used to barcode individual samples; usually CD45 for 'live-cell' barcoding. A 'barcode_node' column must also be present -- the result of \code{FlowSOM::SOM$mapping}.
#' @param barcode.dims Unless defined, \code{barcode.dims} will include any 'CD45' columns -- those used for 'live-cell' barcoding and as input for \code{FlowSOM::SOM}
#' @param m Numeric. Default = 3. Argument passed to \code{utils::combn(x,m)}; used to define the barcode 'scheme': x-choose-m
#' @param delta Numeric. Default = 0.2. Defines the delta/minimum distance between lowest positive barcode marker and highest negative barcode marker; used when assigning barcodes to 'ambiguous' nodes -- those that do not have exactly 3 positive markers by the 'valley' metric.
#' @param ... Argument passed to \code{SOMnambulate:::get.valley(...)}; specifically: \code{no.valley.override}
#'
#' @return Numeric vector of barcode assignments.
#' @export
#'
#'
barcode.assignment<-function(dat,barcode.dims=NULL,m=3,delta=0.2,...){
  ##
  if(!data.table::is.data.table(dat)) stop("Need a data.table")
  if(!'barcode_node' %in% names(dat)) stop("Need a 'barcode_node' column")
  if(is.null(barcode.dims)) barcode.dims<-grep("CD45_",names(dat),value = T)
  l<-length(barcode.dims)
  ##
  means<-dat[,lapply(.SD,mean),.SDcols = barcode.dims,keyby=barcode_node]
  valleys<-dat[,lapply(.SD,get.valley,...),.SDcols = barcode.dims]
  barcode.assignment<-apply(means[,barcode.dims,with = F],1,function(x,v=valleys,key=utils::combn(l,m)){
    s<-sum(x>v)
    if(s>=4){
      0
    }else if(s==3){
      which(apply(key,2,function(i) all(i==which(x>v))))
    }else if(s<=2){
      if((sort(x)[[l-2]]-sort(x)[[l-3]])>delta){
        which(apply(key,2,function(i) all(i==sort(order(x,decreasing = T)[1:3]))))
      }else{
        0
      }
    }
  })
}
#' @title Plot \code{fsom$codes} valleys
#'
#' @param fsom Object; the result of \code{FlowSOM::SOM(...)}
#' @param out.name Character string for use in plot title; usually the 'batch' name
#' @param plot.codes Numeric vector; if defined, will indicate expression level of specific codes
#'
#' @return \code{ggplot2} object
#'
#'
#'
valleys.density.rug<-function(fsom,out.name,plot.codes=NULL){
  ##
  valleys<-data.frame(variable=colnames(fsom$codes),
                      valley=apply(fsom$codes,2,get.valley),
                      # valley.sd=apply(fsom$codes,2,get.valley,sd.adjust=2),
                      row.names = NULL
  )
  dat.codes<-data.table::as.data.table(cbind(fsom$codes,node=seq(nrow(fsom$codes))))
  ##
  data.table::setorder(dat.codes,node)
  dat.codes<-data.table::melt(dat.codes,id.vars=c('node'))
  ##
  valleys.plot<-ggplot2::ggplot() +
    ggplot2::geom_rect(data=valleys,ggplot2::aes(xmin = -Inf, xmax = valley, ymin = -Inf, ymax = Inf),
                       alpha = 1/15,fill = "red") +
    ggplot2::geom_rect(data=valleys,ggplot2::aes(xmin = valley, xmax = Inf, ymin = -Inf, ymax = Inf),
                       alpha = 1/15,fill = "green") +
    ggplot2::geom_density(data=dat.codes,ggplot2::aes(x = value),fill="gray",alpha=.2) +
    ggplot2::geom_rug(data=dat.codes,ggplot2::aes(x = value)) +
    ggplot2::geom_vline(data = valleys,ggplot2::aes(xintercept = valley),color="red",linewidth=1,linetype='dashed') +
    # ggplot2::geom_vline(data = valleys,ggplot2::aes(xintercept = valley.sd),color="purple",linewidth=1,linetype='dashed') +
    ggplot2::geom_text(data = valleys, ggplot2::aes(x = valley, label = round(valley,3)), y = Inf,
                       hjust = -0.25, vjust = 1.5) +
    ggplot2::facet_wrap(~variable) +
    ggplot2::theme(panel.background = ggplot2::element_blank()) +
    ggplot2::xlab("Expression (Scaled)") +
    ggplot2::labs(
      subtitle = out.name,
      title = "Density Distribution of FlowSOM Nodes Per CD45 Expression",
      caption = paste("Density distribution of FlowSOM nodes (n=100) per unique CD45;",
                      "rug-ticks indicate expression value of individual nodes;",
                      "vertical, dashed red-line indicates derived 'valley' value for assignment of a 0 or 1 'barcode key-value'.",
                      sep="\n")
    )
  if(!is.null(plot.codes)){
    valleys.plot <- valleys.plot +
      ggplot2::geom_vline(data = dat.codes[node %in% plot.codes],ggplot2::aes(xintercept= value),
                          color="blue",linewidth=0.25) +
      ggplot2::labs(subtitle = paste(valleys.plot$labels$subtitle,paste("Node(s):",paste(plot.codes,collapse = ", ")),sep="\n"))
  }
  return(valleys.plot)
}
#' @title Plot \code{fsom$codes} valleys with barcode-specific/code-specific indication of expression levels
#'
#' @param fsom Object; the result of \code{FlowSOM::SOM(...)}
#' @param out.name Character string for use in plot title; usually the 'batch' name
#' @param barcode.ids Named character vector; names should equal that of all unique barcode values, including zero; elements should be 'sample.ids'
#'
#' @return a printable list containing barcode-specific \code{ggplot2} objects
#' @export
#'
#'
valleys.density.rug.barcodes<-function(fsom,out.name,barcode.ids=NULL){
  ##
  barcode<-barcode.id<-node<-valley<-value<-NULL
  ##
  barcode.dims<-colnames(fsom$codes)
  valleys<-data.frame(variable=barcode.dims,
                      valley=apply(fsom$codes,2,get.valley),
                      # valley.sd=apply(fsom$codes,2,get.valley,sd.adjust=2),
                      row.names = NULL
  )
  dat.codes<-data.table::as.data.table(cbind(fsom$codes,node=seq(nrow(fsom$codes)),barcode=fsom$barcodes))
  if(!is.null(barcode.ids)){
    dat.codes[,barcode.id:=barcode.ids[as.character(barcode)]]
  }
  ##
  data.table::setorder(dat.codes,barcode)
  dat.codes<-data.table::melt(dat.codes,measure.vars=barcode.dims)
  ##
  plotlist<-lapply(split(dat.codes[barcode!=0],by='barcode'),function(i){
    barcode.num<-i[,unique(barcode)]
    if(!is.null(barcode.ids)){
      barcode.id<-i[,unique(barcode.id)]
    }
    barcode.dims.pos<-barcode.dims[utils::combn(length(barcode.dims),3)[,barcode.num]]
    barcode.alias<- paste0(stringr::str_extract(barcode.dims.pos,"[0-9]{3}"),collapse = " : ")
    # codes.reassigned<-i[,unique(node)][sapply(i[,unique(node)],function(x){
    #   any(i[node==x&variable %in% barcode.dims.pos,value]<=thresholds[barcode.dims.pos])
    # })]
    codes.reassigned<-NULL
    ##
    ggplot2::ggplot() +
      ggplot2::geom_rect(data=valleys,ggplot2::aes(xmin = -Inf, xmax = valley, ymin = -Inf, ymax = Inf),
                         alpha = 1/15,fill = "red") +
      ggplot2::geom_rect(data=valleys,ggplot2::aes(xmin = valley, xmax = Inf, ymin = -Inf, ymax = Inf),
                         alpha = 1/15,fill = "green") +
      ggplot2::geom_density(data=dat.codes,ggplot2::aes(x = value),fill="gray",alpha=.5) +
      ggplot2::geom_rug(data=dat.codes,ggplot2::aes(x = value)) +
      ggplot2::geom_vline(data = valleys,ggplot2::aes(xintercept = valley),color="red",linewidth=0.25,linetype='dashed') +
      ggplot2::geom_vline(data = i,ggplot2::aes(xintercept= value),color="blue",linewidth=0.25) +
      ggplot2::facet_wrap(~variable) +
      ggplot2::theme(panel.background = ggplot2::element_blank()) +
      ggplot2::xlab("Expression (Scaled)") +
      ggplot2::labs(
        title = paste(paste("Barcode#",barcode.num,sep="  "),"Density Distribution of FlowSOM Nodes Per CD45 Expression",sep=" ; "),
        subtitle = paste(
          paste0("Batch: ",out.name),
          if(!is.null(barcode.ids)){
            paste(paste("Barcode Alias:",barcode.alias),
                  paste("Barcode ID:",barcode.id),sep = "    ")
          }else{
            paste("Barcode Alias:",barcode.alias)
          },
          if(length(codes.reassigned)>0){
            paste(paste("Barcode-specific Node(s) ID:",paste0(i[,unique(node)],collapse = ",")),
                  paste("Reassigned Node(s) ID:",paste0(codes.reassigned,collapse = ",")),sep = "    ")
          }else{
            paste("Barcode-specific Node(s) ID:",paste0(i[,unique(node)],collapse = ", "))
          },
          sep ="\n"),
        caption = paste("Density distribution of FlowSOM nodes (n=100) per unique CD45;",
                        "rug-ticks indicate expression value of individual nodes;",
                        "vertical, dashed red-line indicates derived threshold value for assignment of a 0 (red-shade) or 1 (green-shade) 'barcode' value;",
                        "vertical, solid blue-line indicates expression value of barcode-specific nodes",
                        sep="\n")
      )
  })
}
#' @title Plot barcode yields
#'
#' @param dat a \code{data.table} with a 'barcode' column
#' @param out.name Character string for use in plot title; usually the 'batch' name
#' @param barcode.ids Named character vector; names should equal that of all unique barcode values, including zero; elements should be some form of unique identifier (usually 'sample.id')
#'
#' @return a printable list of \code{ggplot2} objects; barcode-specific yields: absolute number and '%'
#' @export
#'
#'
barcode.yield.plot<-function(dat,out.name,barcode.ids=NULL){
  ##
  if(!data.table::is.data.table(dat)&'barcode' %in% names(dat)){
    stop("Need a data.table with a 'barcode' column")
  }
  n<-rbind(
    cbind(dat[,.N,keyby=.(barcode)],type='Barcode 0 Included'),
    cbind(dat[barcode!=0,.N,keyby=.(barcode)],type='Barcode 0 Excluded')
  )
  n[,type := factor(type,levels=c('Barcode 0 Included','Barcode 0 Excluded'))]
  n[,y := N/sum(N)*100,by=.(type)]
  if(!is.null(barcode.ids)){
    n[,barcode.ids:=barcode.ids[as.character(barcode)]]
  }
  n[,barcode:=factor(barcode)]
  ##
  plot.list<-sapply(c('yield','N'),function(i){
    p<-ggplot2::ggplot(n,ggplot2::aes(x=!!ggplot2::sym(ifelse(is.null(barcode.ids),'barcode','barcode.ids')),y=!!ggplot2::sym(ifelse(i=='yield','y','N')))) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~type,scales = "free_y") +
      ggplot2::geom_text(data=n[,.(total=sum(N)),by=.(type)],ggplot2::aes(x=Inf,y=Inf,hjust=1.15,vjust=1.15,label=paste("Total events:",total))) +
      ggplot2::xlab("Barcode #") +
      ggplot2::ylab(ifelse(i=='yield',"Yield (% of total)","# of Events")) +
      ggplot2::labs(
        title=paste0("Batch: ",out.name)
        #caption=paste("Total events:",n[,sum(N)])
      )
    if(!is.null(barcode.ids)){
      p<-p+ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=1,vjust=0.5))
    }
    return(p)
  },simplify = F)
}
#' @title Plot bivariate pairs showing sample-specific de-barcoding results
#'
#' @param dat a \code{data.table} with a 'barcode' column
#' @param out.name Character string for use in plot title; usually the 'batch' name
#' @param barcode.ids Named character vector; names should equal that of all unique barcode values, including zero; elements should be some form of unique identifier (usually 'sample.id')
#' @param bin.number Numeric; argument passed to \code{ggplot2::geom_hex(...,bins=bin.number)}
#'
#' @return a list of \code{gridExtra} objects containing arranged \code{ggplot2} objects
#' @export
#'
#'
barcode.sample.plots<-function(dat,out.name,barcode.ids,bin.number=100){
  ##
  if (!"barcode" %in% names(dat)) {
    stop("Need a 'barcode' column")
  }
  if(!all(data.table::key(dat) == c('barcode','barcode_node'))) data.table::setorder(dat,barcode,barcode_node)
  ##
  barcode.dims <- grep("CD45_", names(dat), value = T)
  barcodes<-dat[,unique(barcode)]
  barcode.key <- utils::combn(length(barcode.dims),3)
  ##
  mdat<-dat[,.N,by=.(barcode)]
  mdat[,yield:=round(N/mdat[,sum(N)]*100,2),by=.(barcode)]
  mdat[,barcode.id:=barcode.ids[as.character(barcode)]]
  ##
  quantile.trim<-sort(unique(dat[,unlist(lapply(.SD,function(x) .I[x>stats::quantile(x,.99999)]),use.names = F),.SDcols=barcode.dims]))
  plot.lims <- dat[-quantile.trim, lapply(.SD, function(x) range(x)),.SDcols=barcode.dims]
  pair.list<-plot.pairs(barcode.dims)
  ##
  dat.bkgd<-dat[-quantile.trim][{set.seed(1337);sample(.N,1E5)}]
  ##
  barcode.plot.list<-lapply(split(dat,by='barcode'),function(dat.bc){
    bc<-dat.bc[,unique(barcode)]
    ##
    plot.list <- lapply(pair.list, function(i, cartesian_lims = plot.lims) {
      p <- ggplot2::ggplot(dat.bkgd, ggplot2::aes(x = !!ggplot2::sym(i[1]),
                                                  y = !!ggplot2::sym(i[2])))
      p <- p + ggplot2::geom_hex(fill = "gray", bins = bin.number)
      p <- p + ggplot2::geom_hex(data = dat.bc,bins = bin.number) +
        viridis::scale_fill_viridis(option = "plasma",limits = c(0, bin.number), oob = scales::squish) +
        ggplot2::theme_classic() + ggplot2::guides(fill = "none")
      p <- p + ggplot2::coord_cartesian(xlim = cartesian_lims[,get(i[1])],
                                        ylim = cartesian_lims[,get(i[2])])
      return(p)
    })
    title.sample<-paste("Sample:",mdat[barcode==bc,barcode.id])
    title.sample.n<-paste("n =",mdat[barcode==bc,N])
    title.sample.yield<-paste0("(",mdat[barcode==bc,yield],"%)")
    title.barcode.number<-paste("Barcode:",bc)
    title.barcode.combination<-paste(stringr::str_extract(barcode.dims[barcode.key[,bc]],"[0-9]{3}"),collapse=" : ")
    title.barcode<-paste(title.barcode.number,title.barcode.combination,sep="     ")
    ##
    plots.arragned<-gridExtra::arrangeGrob(grobs=plot.list,
                                           nrow=2,
                                           ncol=2,
                                           top=paste(out.name,title.sample,title.barcode,sep="\n"),
                                           bottom=paste(title.sample.n,title.sample.yield,sep=" ; ")
    )
    return(plots.arragned)
  })
}
#' @title Plot bivariate pairs showing convoluted de-barcoding results
#'
#' @param dat0 a \code{data.table} with 'barcode' and 'barcode_node' columns; a subset of data contaning only barcode 'zero' values \code{dat[barcode==0]}
#' @param out.name Character string for use in plot title; usually the 'batch' name
#' @param bin.number Numeric; argument passed to \code{ggplot2::geom_hex(...,bins=bin.number)}
#'
#' @return a list of \code{gridExtra} objects containing arranged \code{ggplot2} objects
#' @export
#'
#'
barcode0.node.plots<-function(dat0,out.name,bin.number=100){
  ##
  if(dat0[, unique(barcode)]!=0) stop("Expect only barcode 'zero' data; did you subset the dat? dat0=dat[barcode==0]")
  ##
  if(!"barcode_node" %in% names(dat0)) stop("Need a 'barcode_node' column")
  ##
  barcode.dims<-grep("CD45_",names(dat0),value=T)
  ##
  mdat<-dat0[,.N,keyby=.(barcode_node)]
  mdat[,yield:=round(N/mdat[,sum(N)]*100,2),by=.(barcode_node)]
  quantile.trim<-sort(unique(dat0[,unlist(lapply(.SD,function(x) .I[x>stats::quantile(x,.999)]),use.names = F),.SDcols=barcode.dims]))
  plot.lims <- dat0[-quantile.trim, lapply(.SD, function(x) range(x)),.SDcols=barcode.dims]
  pair.list<-plot.pairs(barcode.dims)
  dat.bkgd<-dat0[-quantile.trim][{set.seed(1337);sample(.N,1E5)}]
  barcode.node.plot.list<-lapply(split(dat0,by="barcode_node",sorted=T),function(dat.bc){
    bc<-dat.bc[,unique(barcode_node)]
    plot.list<-lapply(pair.list,function(i,cartesian_lims=plot.lims){
      p <- ggplot2::ggplot(dat.bkgd, ggplot2::aes(x = !!ggplot2::sym(i[1]),
                                                  y = !!ggplot2::sym(i[2])))
      p <- p + ggplot2::geom_hex(fill = "gray", bins = bin.number)
      p <- p + ggplot2::geom_hex(data = dat.bc, bins = bin.number) +
        viridis::scale_fill_viridis(option = "plasma",
                                    limits = c(0, bin.number), oob = scales::squish) +
        ggplot2::theme_classic() + ggplot2::guides(fill = "none")
      p <- p + ggplot2::coord_cartesian(xlim = cartesian_lims[,
                                                              get(i[1])], ylim = cartesian_lims[, get(i[2])])
      return(p)
    })
    title.sample<-'Sample: Convoluted'
    title.sample.n<-paste("n =",mdat[barcode_node==bc,N])
    title.sample.yield<-paste0("(",mdat[barcode_node==bc,yield],"%)")
    title.barcode<-paste("Barcode Node:",bc)
    plots.arragned<-gridExtra::arrangeGrob(grobs=plot.list,
                                           nrow=2,
                                           ncol=2,
                                           top=paste(out.name,title.sample,title.barcode,sep="\n"),
                                           bottom=paste(title.sample.n,title.sample.yield,sep=" ; ")
    )
    return(plots.arragned)
  })
}
