#' @title Density distributions of selected columns from a `data.table`
#' @description
#' Density distributions are pre-calculated using `data.table` methods; usually passed to \link[ggplot2]{ggplot2} for plotting.
#'
#' @param dt `data.table` containing numeric columns and at least one factored column.
#' @param sd.cols Character vector; used to define `data.table`'s \link[data.table:data.table]{.SDcols} argument; if left undefined, will calculate density distributions for all numeric columns.
#' @param by.cols Character vector; used to define `data.table`'s \link[data.table:data.table]{by} argument; if left undefined, will group by all non-numeric columns.
#' @param trim.quantile Numeric vector of length 2; used to define lower/upper quantile \link[stats:quantile]{probabilities}; trims 'extreme events' before calculating density distributions.
#'
#' @return A `data.table` with pre-calculated 'density.x' and 'density.y' columns; used for plotting.
#' @examples
#' dt<-SOMnambulate:::prepared.examples('dt')
#'
#' #a few select columns (numeric) to asinh transform and plot
#' cols.transform<-c('CD3_170Er','CD4_145Nd','CD8a_162Dy','TNFa_152Sm')
#' for(j in cols.transform){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/5))}
#'
#' #subset 'dt' to calculate density distributions by:
#' #'batch', 'stim.condition', and 'aliquot.seq'
#' dt.dens<-dt.density(
#' dt[,.SD,.SDcols=c(cols.transform,'batch','stim.condition','aliquot.seq')],
#' trim.quantile=c(0.001,0.999))
#' dt.dens[]
#'
#' #make a plot list
#' pl<-sapply(dt.dens[,levels(stim.condition)],function(i){
#' ggplot2::ggplot(dt.dens[stim.condition==i],
#' ggplot2::aes(density.x,density.y,color=aliquot.seq)) +
#' ggplot2::geom_line() +
#' ggplot2::facet_wrap(~variable,scales='free') +
#' ggplot2::labs(
#' title=dt.dens[,as.character(unique(batch))],
#' subtitle=i)
#' },simplify=FALSE)
#'
#' #plot results
#' pl
#'
#' @export
#'
dt.density<-function(dt,sd.cols=NULL,by.cols=NULL,trim.quantile=NULL){
  if(is.null(sd.cols)){sd.cols<-dt[,names(.SD),.SDcols=is.numeric]}
  if(is.null(by.cols)){by.cols<-dt[,names(.SD),.SDcols=!is.numeric]}
  dens.list<-dt[,lapply(.SD,function(x){
    list(stats::density(
      if(!is.null(trim.quantile)){
        q.cut<-stats::quantile(x,trim.quantile)
        x<-x[x>q.cut[1]&x<q.cut[2]]
      }else{
        x
      })[c('x','y')])
  }),
  by=by.cols,.SDcols=sd.cols]
  dt.dens<-data.table::rbindlist(sapply(sd.cols,function(i){
    cbind(dens.list[,stats::setNames(sapply(get(i),'[',c('x','y')),c('density.x','density.y')),by=by.cols],variable=i)
  },simplify = F))
  dt.dens[,'variable' := factor(get('variable'),levels=unique(get('variable')))]
  return(dt.dens)
}

get.dt.col.classes<-function(dt){
  col.classes<-unlist(sapply(dt,class))
  col.classes<-sapply(unique(col.classes),function(i){
    names(col.classes[col.classes %in% i])
  },simplify = F)
  return(col.classes)
}

subsample_dt_sample.id<-function(dat,subsample.val){
  dat[dat[,list(I=if(.N>subsample.val){.I[sample(.N,subsample.val)]}else{.I}),by=get('sample.id')][['I']]]
}
