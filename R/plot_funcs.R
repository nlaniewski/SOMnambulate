density.dt<-function(dt,sd.cols,by.cols){
  dens.list<-dt[,lapply(.SD,function(x) list(stats::density(x[x>0.1&x<stats::quantile(x[x>0.1],.99)])[c('x','y')])),by=by.cols,.SDcols = sd.cols]
  dt.dens<-data.table::rbindlist(sapply(sd.cols,function(i){
    cbind(dens.list[,stats::setNames(sapply(get(i),'[',c('x','y')),c('density.x','density.y')),by=by.cols],variable=i)
  },simplify = F))
  dt.dens[,'variable' := factor(get('variable'),levels=unique(get('variable')))]
  return(dt.dens)
}
#' @title Plot overlaid (by factor) density distributions of selected columns from a \code{data.table}
#' @description
#' Density distributions are pre-calculated using \code{data.table} methods then passed to \code{ggplot2} for plotting.
#'
#' @param dt \code{data.table} containing numeric columns and at least one factored column.
#' @param sd.cols Character vector; used to define \code{data.table}'s \code{.SDcols} argument.
#' @param by.cols Character vector; used to define \code{data.table}'s \code{by} argument; usually just a single factored column for use with \code{ggplot2::facet_wrap()}.
#'
#' @return A \code{ggplot2} object; if no LHS assignment, directly prints.
#' @export
#'
gg.func.density.overlay<-function(dt,sd.cols,by.cols){
  dt.dens<-density.dt(dt,sd.cols,by.cols)
  ggplot2::ggplot(dt.dens,ggplot2::aes(!!as.name('density.x'),!!as.name('density.y'),color=get(names(Filter(is.factor,dt.dens[,!'variable']))))) +
    ggplot2::geom_line(linewidth=0.25) +
    ggplot2::facet_wrap(~variable,scales="free") +
    ggplot2::labs(color=names(Filter(is.factor,dt.dens[,!'variable'])))
}
