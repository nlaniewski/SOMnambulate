barcode.plots <- function(dat,barcode.dims=NULL,bin.number=100){
  if(!'barcode' %in% names(dat)){
    stop("Need a 'barcode' column")
  }
  if(!'pool.id' %in% names(dat)){
    stop("Need a 'pool.id' column")
  }
  if(is.null(barcode.dims)){
    barcode.dims <- grep("CD45_",names(dat),value = T)
  }
  ##
  batch <- dat[,unique(stringr::str_extract(unique(get('pool.id')),"[0-9]{3}_[0-9]{8}"))]
  condition.name <- dat[,ifelse(all(grepl("SEB",unique(get('pool.id')))),"SEB","UNSTIM")]
  batch.name <- paste(batch,condition.name,sep="_")
  if(length(batch.name==1)){
    batch.name <- paste("ECHO",batch.name,sep = "_")
  }else{
    stop("Non-unique 'batch.name'")
  }
  ##
  barcodes <- dat[,sort(unique(get('barcode')))]
  barcode.key <- utils::combn(length(barcode.dims),3)
  ##
  pair.1 <- barcode.dims[seq(1,length(barcode.dims),2)]
  pair.2 <- barcode.dims[seq(2,length(barcode.dims),2)]
  if(length(pair.1)!=length(pair.2)){
    pair.2[length(pair.1)] <- pair.2[length(pair.1)-1]
  }
  pair.list<-mapply(c,pair.1,pair.2,SIMPLIFY = F,USE.NAMES = F)
  ##
  lims <- dat[,lapply(.SD,function(x) c(min(x),max(x))),.SDcols = barcode.dims]
  ##
  barcode.plot.list <- lapply(barcodes,function(bc){
    plot.list <- lapply(pair.list,function(i,cartesian_lims=lims){
      p <- ggplot2::ggplot(dat[sample(.N,1E5),],ggplot2::aes(
        x=!!ggplot2::sym(i[1]),
        y=!!ggplot2::sym(i[2])
      ))
      p <- p + ggplot2::geom_hex(fill = 'gray', bins = bin.number)
      p <- p + ggplot2::geom_hex(data=dat[barcode==bc],bins=bin.number) +
        viridis::scale_fill_viridis(option = "plasma", limits = c(0,bin.number), oob = scales::squish) +
        ggplot2::theme_classic() +
        ggplot2::guides(fill='none')
      p <- p +
        ggplot2::coord_cartesian(
          xlim = cartesian_lims[,get(i[1])],
          ylim = cartesian_lims[,get(i[2])]
        )
      return(p)
    })
    title.sample <- paste("Sample:",dat[barcode==bc,unique(sample)])
    title.sample.n <- paste("n =",dat[barcode==bc,.N])
    title.sample.yield<-paste0("(",round(dat[barcode==bc,.N]/dat[,.N]*100,2),"%)")
    title.barcode.number<-paste("Barcode:",bc)
    title.barcode.combination<-paste(stringr::str_extract(barcode.dims[barcode.key[,bc]],"[0-9]{3}"),collapse = " : ")
    title.barcode<-paste(title.barcode.number,title.barcode.combination,sep = "     ")
    plots.arragned <- gridExtra::arrangeGrob(
      grobs=plot.list,
      nrow=2,
      ncol=2,
      top=paste(batch.name,
                title.sample,
                title.barcode,
                sep = "\n"
      ),
      bottom = paste(title.sample.n,title.sample.yield,sep = " ; ")
    )
    return(plots.arragned)
  })
  ##
}
