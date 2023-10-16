#' @title Generate a barcode key
#'
#' @param barcode.dims Character vector
#' @param m argument passed to \code{combn(x,m)}
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
