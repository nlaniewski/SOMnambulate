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

#' @title Get .N counts and frequencies (of 100%) by grouping variable(s)
#' @description
#' Get \link[data.table]{.N} counts and frequencies (of 100%) \link[data.table:data.table]{by} grouping variable(s). Following \link{map.som.data}, the `node` and `cluster` columns in `dt` can be grouped by (usually) `sample.id` to calculate counts/frequencies of; for use in plots/downstream analyses.
#'
#' @param dt.subset A subset `dt` as returned from \link{map.som.data}; use \link[data.table]{data.table} syntax: `dt[i=NULL,.N,by=...]`
#'
#' @return a 'keyed' `data.table` with `N` and `N.freq` values, grouped `by` columns included in `dt.subset`
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples(example.type='dt')
#'
#' pbmc.markers<-c("CD3","CD4","CD8a","CD14","CD19","CD56","TCRgd")
#' pbmc.dims<-grep(paste0(pbmc.markers,"_",collapse="|"),names(dt),value = TRUE)
#' for(j in pbmc.dims){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#'
#' fsom<-SOMnambulate:::prepared.examples(example.type='fsom')
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=10,umap.codes = TRUE)
#'
#' #do not assign to environment; the function updates 'dt' by reference
#' map.som.data(fsom,dt)
#'
#' #cluster counts
#' dt.N.counts(dt[,.N,keyby=.(cluster_pbmc)])[]
#' dt.N.counts(dt[,.N,keyby=.(sample.id,cluster_pbmc)])[]
#'
#' #node counts
#' dt.N.counts(dt[,.N,keyby=.(node_pbmc)])[]
#' dt.N.counts(dt[,.N,keyby=.(sample.id,node_pbmc)])[]
#'
#' #calculate 'global' node 'sizes'
#' node.counts<-dt.N.counts(dt[,.N,keyby=.(sample.id,stim.condition,node_pbmc)])
#' node.sizes<-node.counts[,list(node.size.median=stats::median(N.freq)),keyby=.(node_pbmc)]
#'
#' #merge and use for plots
#' ggplot2::ggplot(data=merge(fsom$codes.dt,node.sizes)) +
#' ggplot2::aes(umap.1,umap.2,size=node.size.median) +
#' ggplot2::theme_void() +
#' ggplot2::theme(panel.grid = ggplot2::element_blank()) +
#' ggplot2::geom_point(ggplot2::aes(color=cluster_pbmc),alpha=.5) +
#' ggplot2::scale_size_area(max_size = 15,guide="none") +
#' ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3)))
#'
#' #melt, merge and use for plots
#' codes.melt<-data.table::melt(
#' data=fsom$codes.dt,
#' id.vars=grep('node|cluster|umap',names(fsom$codes.dt)))
#'
#' codes.melt[,
#'   value.scaled:= (value-min(value))/(max(value)-min(value)),
#'   by=variable]
#'
#' ggplot2::ggplot(data=merge(codes.melt,node.sizes)) +
#' ggplot2::aes(umap.1,umap.2,color=value.scaled,size=node.size.median) +
#' ggplot2::geom_point(alpha=.75) +
#' ggplot2::scale_size_area(max_size = 10,guide="none") +
#' viridis::scale_color_viridis(option='plasma') +
#' ggplot2::theme_dark() +
#' ggplot2::theme(panel.grid = ggplot2::element_blank()) +
#' ggplot2::facet_wrap(~variable)
#'
dt.N.counts<-function(dt.subset){
  .by<-names(dt.subset)[!names(dt.subset) %in% "N"]
  if(length(.by)>1){.by<-.by[1]}else{.by<-NULL}
  dt.subset[,total:=sum(N),by=.by
  ][,N.freq:=(N/total)*100
  ][,total:=NULL
  ]
}
