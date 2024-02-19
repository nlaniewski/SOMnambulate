##file.paths
extdata<-system.file("extdata",package="SOMnambulate")
fcs.files <- list.files(extdata,full.names=TRUE,pattern=".fcs")

##prepare 'fcs.files.dt_ECHO' for use in examples
fcs.files.dt<-get.fcs.file.dt(fcs.files)
new.cols<-c('study.name','batch.seq','batch.date')
fcs.files.dt[,(new.cols) :=
data.table::tstrsplit(basename(dirname(FILENAME)),"_",type.convert=list(as.factor=1:3))]
fcs.files.dt[,stim.condition := factor(stringr::str_extract(FIL,"SEB|UNSTIM"))]
fcs.files.dt[,aliquot.seq := factor(seq(.N)),by=.(batch.date,stim.condition)]
fcs.files.dt[,batch := factor(paste(study.name,batch.seq,batch.date,sep="_"))]
fcs.files.dt[,sample.id := paste(batch,stim.condition,aliquot.seq,sep="_")]

fcs.files.dt_ECHO<-fcs.files.dt

##prepare 'ca_ECHO' for use in examples
ca<-get.fcs.channel.alias(fcs.files)

drop.metals<-c("Y","Pd","Sn","I","La","Lu","Os","Ir","Pt","Pb")
drop.metals<-paste0(drop.metals,'$',collapse="|")
drop.pattern<-paste(drop.metals,'background','noise',sep="|")

ca<-ca[grep(drop.pattern,alias,ignore.case=TRUE,invert=TRUE)]

ca[grep("beads",alias,ignore.case = TRUE),alias:=sub("Norm_beads","beads",alias)]
ca[grep("viability",alias,ignore.case = TRUE),alias:="194Pt_viability"]#drop '_cisplatin'

ca[stringr::str_detect(alias,"[A-Z]{1}[a-z]{1}_"),alias := sub('(\\w+)_(\\w+)', '\\2_\\1', alias)]
alias.order<-c('Time',sort(grep('event|center|offset|width|residual',ca$alias,value = T,ignore.case = T)))
alias.order<-c(alias.order,ca[!alias %in% alias.order,stringr::str_sort(alias,numeric = T)])
data.table::setorder(ca[, 'ord' := match(alias,alias.order)],'ord')[,'ord' := NULL]

ca_ECHO<-ca

usethis::use_data(fcs.files.dt_ECHO,ca_ECHO,overwrite = TRUE,internal = TRUE)
