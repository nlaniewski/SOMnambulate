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

barcode.assignment.codes<-function(codes,k=3){
  n<-ncol(codes)
  barcode.combinations <- t(utils::combn(n,k))
  barcode.key <- matrix(data = 0,
                        nrow = nrow(barcode.combinations),
                        ncol = n,
                        dimnames = list(NULL, colnames(codes))
  )
  for(i in seq(nrow(barcode.key))){
    barcode.key[i, barcode.combinations[i, ]] <- 1
  }
  node.thresholds <- apply(codes, 2, function(x){
    q.vals <- stats::quantile(x,probs=c(.001,.999))
    d <- stats::density(x[x>q.vals[1]&x<q.vals[2]])
    valley.max.x = max(d$x[which(diff(sign(diff(d$y)))==2)])
    return(valley.max.x)
  })
  node.key <- sapply(names(node.thresholds), function(x){
    key <- ifelse(codes[, x] > node.thresholds[x], 1, 0)
  })
  node.key.t <- t(node.key)
  barcode.key.t <- t(barcode.key)
  barcode.assignment <- apply(node.key.t, 2, function(i) {
    if(max(colSums(i == barcode.key.t)) == n){
      barcode <- which(colSums(i == barcode.key.t) == n)
    }else{
      barcode <- 0
    }
  })
  return(barcode.assignment)
}
