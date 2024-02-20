?fcs.to.dt.parallel
fcs.files.dt<-prepared.examples('fcs.files.dt')

#from the 'get.fcs.channel.alias' example:
ca<-prepared.examples('channel.alias')

dt<-fcs.to.dt.parallel(
  fcs.file.dt=fcs.files.dt[,.(f.path,sample.id,batch,stim.condition,aliquot.seq)],
  channel_alias=ca,use.alias.pattern=TRUE,alias.order=TRUE)
metals<-grep("\\w+_[0-9]{3}",names(dt),value = T)

for(j in metals){data.table::set(dt,j=j,value = asinh(dt[[j]]/10))}

ggplot2::ggplot(
  dt.density(dt[,.SD,.SDcols = c(metals,'stim.condition')],trim.quantile = c(0.01,0.99)),
  ggplot2::aes(density.x,density.y,color=stim.condition)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~variable,scales="free")

ggplot2::ggplot(
  dt.density(dt[stim.condition=="UNSTIM",.SD,.SDcols = c(metals,'aliquot.seq')],trim.quantile = c(0.01,0.99)),
  ggplot2::aes(density.x,density.y,color=aliquot.seq)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~variable,scales="free")

library(ggplot2)
ggplot(dt,aes(Time,beads_140Ce)) +
  geom_hex(bins=100) +
  geom_hline(data=dt[,.(y.int=break.detection.histogram.counts(beads_140Ce)),by=.(stim.condition)],aes(yintercept = y.int)) +
  facet_wrap(~stim.condition,scales="free")

ggplot(dt,aes(Time,viability_194Pt)) +
  geom_hex(bins=100) +
  geom_hline(data=dt[,.(y.int=break.detection.histogram.counts(viability_194Pt)),by=.(stim.condition)],aes(yintercept = y.int)) +
  facet_wrap(~stim.condition,scales="free")
