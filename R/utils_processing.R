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
##
drop_extreme_events_fluors<-function(x,count.cut=3,quantile.cut=0.05){
  b <- seq(min(x), max(x), length.out = 301)
  h <- graphics::hist(x, breaks = b, plot = F)
  count.cut.val<-floor(stats::quantile(h$counts[h$counts>count.cut],quantile.cut))
  count.range<-range(which(h$counts>count.cut.val))
  break.range<-range(h$breaks[count.range[1]:count.range[2]])
  ##
  drop.index<-c(which(x<break.range[1]),which(x>break.range[2]))
  return(drop.index)
}
##
drop_extreme_events_time<-function(time.vector,plot.histogram=F,sample.name){
  b <- seq(min(time.vector), max(time.vector), length.out = 301)
  h <- graphics::hist(time.vector, breaks = b, plot = F)
  if(any(h$counts<(mean(h$counts)*.05))){
    counts.cut<-min(which(h$counts<(mean(h$counts)*.05)))-1
    hc<-h$counts[1:counts.cut]
  }else{
    hc<-h$counts
  }
  for (i in 5:1) {
    if (sign(stats::median(hc) - stats::sd(hc) *
             i) == (-1)) {
      next
    }
    else {
      drop.bins.i <- which(h$counts < (stats::median(hc) -
                                         stats::sd(hc) * i))
      break
    }
  }
  ##
  if(length(drop.bins.i) == 0){
    drop.bin<-length(h$counts)
    time.break <- max(time.vector)
  }else{
    drop.bin<-min(drop.bins.i)
    time.break <- h$breaks[drop.bin]
  }
  if(abs(which.max(h$counts)-drop.bin)<5){
    time.break<-h$breaks[max(which(sign(diff(h$counts[1:which.max(h$counts)]))==(-1)))]
  }
  ##
  if(plot.histogram){
    plot(h,main=sample.name)
    graphics::abline(v=time.break,col="red",lty='dashed',lwd=2)
  }else{
    return(time.break)
  }
}
##
