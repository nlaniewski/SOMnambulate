#' @title Density distributions of selected columns from a \code{data.table}
#' @description
#' Density distributions are pre-calculated using \code{data.table} methods; usally passed to \code{ggplot2} for plotting.
#'
#' @param dt \code{data.table} containing numeric columns and at least one factored column.
#' @param sd.cols Character vector; used to define \code{data.table}'s \code{.SDcols} argument; if left undefined, will calculate density distributions for all numeric columns.
#' @param by.cols Character vector; used to define \code{data.table}'s \code{by} argument; if left undefined, will group by all non-numeric columns.
#' @param trim.quantile Numeric vector of length 2; used to define \code{stats::quantile}'s \code{probs} argument; trims 'extreme events' before calculating density distributions.
#'
#' @return A \code{data.table} with pre-calculated 'density.x' and 'density.y' columns; used for plotting.
#'
density.dt<-function(dt,sd.cols=NULL,by.cols=NULL,trim.quantile=NULL){
  if(is.null(sd.cols)){sd.cols<-dat[,names(.SD),.SDcols=is.numeric]}
  if(is.null(by.cols)){by.cols<-dat[,names(.SD),.SDcols=!is.numeric]}
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
# subsample_dt_sample.id<-function(dat,subsample.val){
#   dat[dat[,list(I=if(.N>subsample.val){.I[sample(.N,subsample.val)]}else{.I}),by=get('sample.id')][['I']]]
# }
