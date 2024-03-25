#' @title Prepared examples
#'
#' @param example.type character string; matches one of: c('fcs.files.dt','channel.alias','dt','fsom')
#'
#' @return (hopefully) something useful for building examples and reducing redundancy
#'
#'
#' @examples
#' SOMnambulate:::prepared.examples(example.type='fcs.files.dt')[]
#' SOMnambulate:::prepared.examples(example.type='channel.alias')[]
#' SOMnambulate:::prepared.examples(example.type='dt')[]
#' str(SOMnambulate:::prepared.examples(example.type='fsom'))
prepared.examples <- function(example.type=c('fcs.files.dt','channel.alias','dt','fsom')){
  #R CMD check; NULL data.table vars
  FILENAME<-stim.condition<-FIL<-aliquot.seq<-batch<-batch.date<-batch.seq<-study.name<-sample.id<-NULL
  #
  extdata<-system.file("extdata",package="SOMnambulate")
  fcs.files <- list.files(extdata,full.names=TRUE,pattern=".fcs")
  ##prepare 'fcs.files.dt' as in 'get.fcs.file.dt' (examples)
  fcs.files.dt_example<-function(){
    message(paste(
      "A data.table 'fcs.files.dt' as prepared by:"
      ,"get.fcs.file.dt()"
      ,"see specific function examples"
      ,sep = "\n"
    ))
    fcs.files.dt<-get.fcs.file.dt(fcs.files)
    new.cols<-c('study.name','batch.seq','batch.date')
    fcs.files.dt[,(new.cols) :=
                   data.table::tstrsplit(basename(dirname(FILENAME)),"_",type.convert=list(as.factor=1:3))]
    fcs.files.dt[,stim.condition := factor(stringr::str_extract(FIL,"SEB|UNSTIM"))]
    fcs.files.dt[,aliquot.seq := factor(seq(.N)),by=.(batch.date,stim.condition)]
    fcs.files.dt[,batch := factor(paste(study.name,batch.seq,batch.date,sep="_"))]
    fcs.files.dt[,sample.id := paste(batch,stim.condition,aliquot.seq,sep="_")]
    #
    return(fcs.files.dt)
  }
  ##prepare 'channel.alias' as in 'get.fcs.channel.alias' (examples)
  channel.alias_example<-function(){
    message(paste(
      "A data.table 'channel_alias' as prepared by:"
      ,"get.fcs.channel.alias()"
      ,"see specific function examples"
      ,sep = "\n"
    ))
    drop.metals<-c("89Y",paste0(c(102,104:106,108,110),"Pd"),"120Sn","127I",
                   paste0(c(138,139),"La"),"176Lu","190Os",paste0(c(191,193),"Ir"),
                   paste0(c(195,198),"Pt"),"208Pb")
    drop.metals<-paste0(drop.metals,'$',collapse="|")
    drop.pattern<-paste(drop.metals,'background','noise',sep="|")
    ca<-get.fcs.channel.alias(fcs.files,drop.pattern.alias = drop.pattern, order.alias = T)
    ca[grep("beads",alias,ignore.case = TRUE),alias:=sub("Norm_beads","beads",alias)]
    ca[grep("viability",alias,ignore.case = TRUE),alias:="194Pt_viability"]#drop '_cisplatin'
    ca[stringr::str_detect(alias,"[A-Z]{1}[a-z]{1}_"),alias := sub('(\\w+)_(\\w+)', '\\2_\\1', alias)]
    alias.order<-c('Time',sort(grep('event|center|offset|width|residual',ca$alias,value = T,ignore.case = T)))
    alias.order<-c(alias.order,ca[!alias %in% alias.order,stringr::str_sort(alias,numeric = T)])
    data.table::setorder(ca[, 'ord' := match(alias,alias.order)],'ord')[,'ord' := NULL]
    #
    return(ca)
  }
  ##prepare 'dt' as in 'fcs.to.dt' (examples)
  dt_example<-function(){
    message(paste(
      "A data.table 'dt' as prepared by:"
      ,"get.fcs.file.dt()"
      ,"get.fcs.channel.alias()"
      ,"fcs.to.dt()"
      ,"see specific function examples"
      ,sep = "\n"
    ))
    suppressMessages(dt<-data.table::rbindlist(
      sapply(
        split(fcs.files.dt_example()[,.(f.path,sample.id,batch,stim.condition,aliquot.seq)],by='f.path'),
        fcs.to.dt,channel_alias=channel.alias_example(),use.alias.pattern=TRUE,alias.order=TRUE,simplify=FALSE)
    ))
  }
  ##prepare 'fsom' as in 'som' (examples)
  som_example<-function(){
    dt<-suppressMessages(dt_example())
    message("'fsom' as prepared by 'som()'")
    pbmc.markers<-c("CD3","CD4","CD8a","CD14","CD19","CD56","TCRgd")
    pbmc.dims<-grep(paste0(pbmc.markers,"_",collapse="|"),names(dt),value = T)
    for(j in pbmc.dims){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
    fsom<-som(dt[,pbmc.dims,with=FALSE],.scale=T,map=FALSE)
    return(fsom)
  }
#
  if(example.type=='fcs.files.dt') return(fcs.files.dt_example())
  if(example.type=='channel.alias') return(channel.alias_example())
  if(example.type=='dt') return(dt_example())
  if(example.type=='fsom') return(som_example())
}
