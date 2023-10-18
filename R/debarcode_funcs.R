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
#' @param sd.adjust Numeric; if defined, will adjust the return value by multiples of density value standard deviations.
#'
#' @return Numeric representing the 'valley' value in a density distribution.
#'
#'
#'
get.valley<-function(x,quantile.trim=T,sd.adjust=NULL){
  if(quantile.trim){
    q.vals <- stats::quantile(x, probs = c(0.001, 0.999))#trim 'extreme' values
    d <- stats::density(x[x > q.vals[1] & x < q.vals[2]])
  }else{
    d <- stats::density(x)
  }
  diff.derivative.1<-sign(diff(d$y))#(sign of) derivative 1 of density 'y' values; monotonic 'slopes'
  diff.derivative.2<-diff(diff.derivative.1)#derivative 2; peaks and valleys
  valley<-which(diff.derivative.2==2)#-2 == peak(s);2 == valley
  if(length(valley)==1){
    if(!is.null(sd.adjust)){
      return(d$x[valley]-stats::sd(d$y)*sd.adjust)#return x value where valley occurs minus 'sd.adjust' multiples of y value standard deviations; 'expression' value
    }else{
      return(d$x[valley])#return x value where valley occurs; 'expression' value
    }
  }else{
    stop(paste(paste(length(valley),"Valleys detected..."),
               paste("Valleys at:",paste(round(d$x[valley],5),collapse = ' ; ')),
               sep='\n')
    )
  }
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
  node<-valley<-value<-NULL
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
  barcode<-total<-type<-y<-NULL
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
