subsample_dt_sample.id<-function(dat,subsample.val){
  dat[dat[,list(I=if(.N>subsample.val){.I[sample(.N,subsample.val)]}else{.I}),by=get('sample.id')][['I']]]
}
