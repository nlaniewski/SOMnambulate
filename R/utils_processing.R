break.detection.histogram.counts<-function(x,zero.peak.knock.down=T){
  if(zero.peak.knock.down){
    x<-x[x>0]
  }
  m<-mean(x)
  x<-x[x>m]
  ##
  b <- seq(min(x), max(x), length.out = 201)
  h <- graphics::hist(x,breaks=b,plot = F)
  ##
  dc<-diff(h$counts)
  med<-stats::median(dc)
  break.i<-which(dc>med)[1]+1
  ##
  return(h$breaks[break.i])
}
