#' @title A cluster-faceted bivariate ggplot with a silhouette sampled from all data points
#'
#' @param dt a `data.table` containing a 'cluster' column
#' @param ... `ggplot2::aes` arguments: x and y (unquoted)
#' @param bins argument as defined by `ggplot2::geom_hex`
#' @param fill.limits Numeric vector of length two used to define the `limits` argument of `ggplot2::scale_fill_gradient(...)`
#' @param silhouette.val Numeric; the number of data points to sample from `dt` to create the silhouette
#'
#' @return a \code{ggplot2} object
#' @export
#'
#'
gg.func.bivariate.cluster.silhouette<-function(dt,...,bins=200,fill.limits=c(0,bins/2),silhouette.val=5E4){
  cluster.col=grep('cluster',names(dt),value = T)
  p<-ggplot2::ggplot(data=NULL,ggplot2::aes(...))
  p<-p+ggplot2::geom_hex(data=dt[sample(.N,silhouette.val),!cluster.col,with=F],fill='gray',bins=bins)
  p<-p+ggplot2::geom_hex(data=dt,bins=bins)
  p<-p+viridis::scale_fill_viridis(option="plasma",limits=fill.limits,oob=scales::squish)
  p<-p+ggplot2::theme_classic()
  p+ggplot2::facet_wrap(~get(cluster.col))
}
