#' @title A bivariate ggplot of hex-binned data points
#'
#' @param dt An object as returned from \link{fcs.to.dt}
#' @param ... \link[ggplot2]{aes} arguments (unquoted); essentially x and y. If a `z` aesthetic is provided, the plot will switch to \link[ggplot2]{stat_summary_hex}, using this defined third variable as a 'color-by'.
#' @param bins \link[ggplot2]{geom_hex} argument; numeric vector giving number of bins in both vertical and horizontal directions.
#' @param limits \link[ggplot2]{continuous_scale} argument.
#'
#' @return a \link[ggplot2]{ggplot} object
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples("dt")
#' metals <- grep("_[0-9]{3}[A-Z]{1}[a-z]{1}",names(dt),value=TRUE)
#' for(j in metals){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#'
#' plotSOMprimary(dt,CD4_145Nd,CD8a_162Dy,bins=100)
#'
#' plotSOMprimary(dt,CD4_145Nd,CD8a_162Dy,z=CD3_170Er,bins=100)
#' plotSOMprimary(dt,CD4_145Nd,CD8a_162Dy,z=CD3_170Er,bins=100,
#'   limits=c(0,dt[,quantile(CD3_170Er,.9)])) +
#'   ggplot2::theme_dark()
plotSOMprimary <- function(dt,...,bins=200,limits=NULL){
  ##
  p<-ggplot2::ggplot(
    data = dt,
    mapping = ggplot2::aes(...)
  )
  if("z" %in% ...names()){
    dot.names=lapply(substitute(list(...))[-1], deparse)
    z<-dot.names[...names()=="z"]
    p <- p +
      ggplot2::stat_summary_hex(bins=bins) +
      viridis::scale_fill_viridis(option="magma",
                                  limits = limits,
                                  oob = scales::squish) +
      ggplot2::labs(fill=z)
  }else{
    p <- p +
      ggplot2::geom_hex(bins=bins) +
      viridis::scale_fill_viridis(option="plasma",
                                  limits = limits,
                                  oob = scales::squish)
  }
  ##
  return(p)
}
#' @title A bivariate ggplot of hex-binned data points -- interactive (shiny)
#'
#' @param dt An object as returned from \link{fcs.to.dt}
#'
#' @return \link[shiny]{shinyApp}
#' @export
#'
#' @examples
#' if (interactive()) {
#' dt<-SOMnambulate:::prepared.examples("dt")
#' metals <- grep("_[0-9]{3}[A-Z]{1}[a-z]{1}",names(dt),value=TRUE)
#' for(j in metals){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#'
#' plotSOMprimary_shiny(dt)
#'}
plotSOMprimary_shiny<-function(dt){
  ##use internal function to get dt column classes
  col.classes<-get.dt.col.classes(dt)
  ##vectors for use as shiny inputs
  choices.numeric<-col.classes$numeric
  choices.factor<-col.classes$factor
  choices.sample<-dt[,sort(unique(sample.id))]
  ##limits for numeric columns; all data
  lims <- dt[,lapply(.SD,range),.SDcols=choices.numeric]
  ##shiny ui inputs
  ui.inputs<-list(
    marker.x=shiny::selectInput(
      inputId = "marker1",
      label = "Marker (x):",
      choices = choices.numeric,
      selected = choices.numeric[1]
    ),
    marker.y=shiny::selectInput(
      inputId = "marker2",
      label = "Marker (y):",
      choices = choices.numeric,
      selected = choices.numeric[2]
    ),
    color.by=shiny::selectInput(
      inputId = "color.by",
      label = "Color:",
      choices = choices.numeric,
      selected = ""
    ),
    color.check=shiny::checkboxInput(
      inputId = "color.check",
      label = "Color by selected variable?",
      value = FALSE
    ),
    sample.id=shiny::selectInput(
      inputId = "sample.id",
      label="Sample:",
      choices=choices.sample,
      selected = choices.sample[1]
    )
  )
  ##shiny ui actionbuttons
  actionbuttons<-lapply(2:4,function(n){
    list(pop=
           shiny::actionButton(
             inputId = paste0("pop",n),
             label = paste0("#",n),
             icon = shiny::icon(name="new-window",lib="glyphicon")
           ),
         clear=
           shiny::actionButton(
             inputId = paste0("clear.pop",n),
             label = paste0("#",n),
             icon = shiny::icon(name="remove",lib="glyphicon")
           )
    )
  })
  ##
  ui <- shinydashboard::dashboardPage(
    header = shinydashboard::dashboardHeader(
      title = "SOMthingShiny",disable = F
    ),
    sidebar = shinydashboard::dashboardSidebar(
      shiny::numericInput(
        inputId = "bins",
        label="Bins:",
        value=200
      ),
      shiny::splitLayout(
        shiny::textInput(
          inputId = "limits.count",
          label = "Limits (count):",
          placeholder = "#,#"
        ),
        shiny::textInput(
          inputId = "limits.color",
          label = "Limits (color):",
          placeholder = "#,#"
        )
      ),
      ui.inputs$marker.x,
      ui.inputs$marker.y,
      ui.inputs$color.by,
      ui.inputs$color.check,
      ui.inputs$sample.id,
      ui.inputs$factors,
      shiny::selectInput(
        inputId = "ggtheme",
        label = "ggtheme:",
        choices = c(
          "gray"
          ,"bw"
          ,"light"
          ,"dark"
          ,"minimal"
          ,"classic"
          ,"void"
        ),
        selected = "dark"
      ),
      lapply(actionbuttons,function(i){shiny::splitLayout(i$pop,i$clear)})
    ),
    body = shinydashboard::dashboardBody(
      shiny::fluidRow(
        lapply(c(".plotSOMprimary","pop2"),function(n){
          shinydashboard::box(
            collapsible = T,
            title = NULL,
            shiny::plotOutput(n)
          )
        })
      ),
      shiny::fluidRow(
        lapply(3:4,function(n){
          shinydashboard::box(
            collapsible = T,
            title = NULL,
            shiny::plotOutput(paste0("pop",n))
          )
        })
      )
    )
  )
  ##
  server = function(input, output) {
    ##function to extract color and count limits (https://stackoverflow.com/users/5836932/k-rohde)
    extract.lims <- function(text) {
      text <- gsub(" ", "", text)
      split <- strsplit(text, ",", fixed = FALSE)[[1]]
      as.numeric(split)
    }
    ##bivariate plot generated by plotSOMprimary
    .plotSOMprimary <- shiny::reactive({
      p <- if(input$color.check){
        plotSOMprimary(
          dt[sample.id %in% input$sample.id],
          x=!!as.name(input$marker1),
          y=!!as.name(input$marker2),
          z=!!as.name(input$color.by),
          bins = input$bins,
          limits = if(!anyNA(extract.lims(input$limits.color)) && length(extract.lims(input$limits.color)) == 2){
            extract.lims(input$limits.color)
          }
        ) +
          ggplot2::labs(fill=input$color.by,subtitle = input$sample.id)
      }else{
        plotSOMprimary(
          dt[sample.id %in% input$sample.id],
          x=!!as.name(input$marker1),
          y=!!as.name(input$marker2),
          bins = input$bins,
          limits = if(!anyNA(extract.lims(input$limits.count)) && length(extract.lims(input$limits.count)) == 2){
            extract.lims(input$limits.count)
          }
        ) +
          ggplot2::labs(subtitle = input$sample.id)
      }
      p <- p +
        eval(parse(text = paste0('ggplot2::theme_',input$ggtheme,'()'))) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank()
        ) +
        #ggplot2::guides(fill='none') +
        if(all(c(input$marker1,input$marker2) %in% names(lims))){
          ggplot2::coord_cartesian(
            xlim=lims[,get(input$marker1)],
            ylim=lims[,get(input$marker2)])
        }else if(input$marker1 %in% names(lims)){
          ggplot2::coord_cartesian(
            xlim=lims[,get(input$marker1)])
        }else if(input$marker2 %in% names(lims)){
          ggplot2::coord_cartesian(
            ylim=lims[,get(input$marker2)])
        }
      return(p)
    })
    ##bivariate plot; output
    output$.plotSOMprimary <- shiny::renderPlot({
      .plotSOMprimary()
    })
    ##there is probably a more efficient way of handling the following code chunks
    ##can probably put these into a list/loop
    ##isolate pop-out #2; the current reactive plot as a new isolated/static plot
    shiny::observeEvent(input$pop2,{
      output$pop2<-shiny::renderPlot({
        shiny::isolate(.plotSOMprimary())
      })
    })
    ##clear pop-out #2
    shiny::observeEvent(input$clear.pop2,{
      output$pop2 <- NULL
    })
    ##isolate pop-out #3; the current reactive plot as a new isolated/static plot
    shiny::observeEvent(input$pop3,{
      output$pop3<-shiny::renderPlot({
        shiny::isolate(.plotSOMprimary())
      })
    })
    ##clear pop-out #3
    shiny::observeEvent(input$clear.pop3,{
      output$pop3 <- NULL
    })
    ##isolate pop-out #4; the current reactive plot as a new isolated/static plot
    shiny::observeEvent(input$pop4,{
      output$pop4<-shiny::renderPlot({
        shiny::isolate(.plotSOMprimary())
      })
    })
    ##clear pop-out #4
    shiny::observeEvent(input$clear.pop4,{
      output$pop4 <- NULL
    })
  }
  ##
  shiny::shinyApp(ui, server)
}
#' @title Density distribution of FlowSOM nodes per unique barcode channel with barcode-specific nodes as rug-ticks
#'
#' @param fsom An object as returned from \link{som.codes.dt}.
#' @param mdat **Optional** data.table of metadata used to associate unique barcodes with unique samples. `mdat` must have a 'barcode' column; all additional columns will be used to construct a plot subtitle.
#'
#' @return a \link[ggplot2]{ggplot} object
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples('dt')
#' dt<-dt[stim.condition=="SEB"]
#' barcode.dims<-grep("CD45_",names(dt),value = TRUE)
#' for(j in barcode.dims){
#' data.table::set(dt,j=j,value=asinh(dt[[j]]/10))
#' }
#' fsom<-som(dt[,barcode.dims,with=FALSE],.scale=FALSE,xdim=9,ydim=9,map=FALSE)
#' fsom<-som.codes.dt(fsom,append.name = 'barcode')
#' SOMnambulate::map.som.data(fsom,dt)
#' fsom$codes.dt[,barcode := key.codes(fsom)$barcode]
#' dt[,barcode := fsom$codes.dt$barcode[node_barcode]]
#'
#' extdata<-system.file("extdata",package="SOMnambulate")
#' mdat <- data.table::fread(
#'   list.files(extdata,full.names = TRUE, pattern="batch_sheet.csv"),
#'   keepLeadingZeros = TRUE
#' )
#' mdat<-mdat[stim.condition=="SEB",.(barcode,sample.id,stim.condition,batch)]
#' plotSOMbarcodes_rug(fsom,mdat)[1:5]
#'
plotSOMbarcodes_rug <- function(fsom,mdat=NULL){
  #get local maxima and minima from fsom$codes
  peaks.valleys<-apply(fsom$codes,2,get.local.maxima.minima)
  #valleys
  valleys<-lapply(peaks.valleys,'[[','lminima')
  #some CD45s may return more than one valley
  #assumption: first element represents the 'true' valley
  valleys<-lapply(valleys,'[[',1)
  #as a data.table -- melted
  valleys<-data.table::melt(
    data.table::as.data.table(valleys),
    measure.vars=names(valleys),
    value.name = 'valley')
  #
  codes.dt.melt<-data.table::melt(fsom$codes.dt,measure.vars=colnames(fsom$codes))
  data.table::setorder(codes.dt.melt,barcode)
  #preserve barcode==0; use all.x=TRUE
  codes.dt.melt<-merge(codes.dt.melt,mdat,all.x = TRUE)
  #
  nodes.n<-codes.dt.melt[,length(unique(node_barcode))]
  barcode.dims<-colnames(fsom$codes)
  #
  plotlist<-lapply(split(codes.dt.melt[barcode!=0],by='barcode'),function(dt.bc){
    #unique barcode number
    bc.num<-dt.bc[,unique(barcode)]
    #unique barcode alias
    bc.alias<-barcode.dims[key.barcode(barcode.dims,3)[bc.num,barcode.dims,with=F]==1]
    bc.alias<-paste0(stringr::str_extract(bc.alias,"_([0-9]+)",group=1), collapse = " : ")
    #unique metadata
    .mdat<-unlist(unique(dt.bc[,names(dt.bc)[!names(dt.bc) %in% c('barcode','node_barcode','variable','value')],with = F]))
    .mdat<-.mdat[sort(names(.mdat))]
    #unique nodes
    nodes<-paste0(dt.bc[,sort(unique(node_barcode))],collapse = ", ")
    #
    ggplot2::ggplot() +
      ggplot2::geom_rect(data=valleys,ggplot2::aes(xmin = -Inf, xmax = valley, ymin = -Inf, ymax = Inf),
                         alpha = 1/15,fill = "red") +
      ggplot2::geom_rect(data=valleys,ggplot2::aes(xmin = valley, xmax = Inf, ymin = -Inf, ymax = Inf),
                         alpha = 1/15,fill = "green") +
      ggplot2::geom_density(data=codes.dt.melt,ggplot2::aes(x = value),fill="gray",alpha=.5) +
      ggplot2::geom_rug(data=codes.dt.melt,ggplot2::aes(x = value)) +
      ggplot2::geom_vline(data = valleys,ggplot2::aes(xintercept = valley),color="red",linewidth=0.25,linetype='dashed') +
      ggplot2::geom_vline(data = dt.bc,ggplot2::aes(xintercept= value),color="blue",linewidth=0.25) +
      ggplot2::facet_wrap(~variable) +
      ggplot2::theme(panel.background = ggplot2::element_blank()) +
      ggplot2::xlab("Expression") +
      ggplot2::labs(
        subtitle = paste(
          paste(paste("Barcode:",bc.num),bc.alias,sep = "     "),
          paste("Barcode-specific Nodes:",nodes),
          "",
          paste(names(.mdat),.mdat,sep = ": ",collapse = "\n"),
          sep="\n"),
        caption = paste(
          paste0("Density distribution of FlowSOM nodes (n=",nodes.n,") per unique CD45;"),
          "rug-ticks indicate expression value of individual nodes;",
          "vertical, dashed red-line indicates derived threshold value for assignment of a 0 (red-shade) or 1 (green-shade) 'barcode' value;",
          "vertical, solid blue-line indicates expression value of barcode-specific nodes",
          sep="\n"))
  })
}
#' @title Bivariate plots of CD45 expression of node-specific barcode 'zero' events
#'
#' @param dt An object as returned from \link{fcs.to.dt} and mapped using \link{map.som.data}.
#' @param fsom An object as returned from \link{som.codes.dt}; additionally, requires a 'barcode' column.
#' @param node.threshold Numeric; to avoid overplotting, if the number of node-specific events exceed the defined threshold value, sampling will take place.
#' @param bins.N Argument as defined by `ggplot2::geom_hex`.
#' @param mdat.cols Character vector; columns in `dt` that will be used to construct a plot subtitle.
#'
#' @return a \code{ggplot2} object
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples('dt')
#' dt<-dt[stim.condition=="SEB"]
#' barcode.dims<-grep("CD45_",names(dt),value = TRUE)
#' for(j in barcode.dims){
#' data.table::set(dt,j=j,value=asinh(dt[[j]]/10))
#' }
#' fsom<-som(dt[,barcode.dims,with=FALSE],.scale=FALSE,xdim=9,ydim=9,map=FALSE)
#' fsom<-som.codes.dt(fsom,append.name = 'barcode')
#' SOMnambulate::map.som.data(fsom,dt)
#' fsom$codes.dt[,barcode := key.codes(fsom)$barcode]
#' dt[,barcode := fsom$codes.dt$barcode[node_barcode]]
#'
#' plotSOMbarcodes_zero(dt,fsom,bins.N=50,mdat.cols=c('batch','stim.condition'))[1:5]
plotSOMbarcodes_zero<-function(dt,fsom,node.threshold=1E4,bins.N=100,mdat.cols=NULL){
  #dimensions used to generate the SOMs/barcodes
  barcode.dims<-colnames(fsom$codes)
  #axis limits
  lims<-dt[,lapply(.SD,range),.SDcols = barcode.dims]
  #plot pairs
  pp<-plot.pairs(barcode.dims)
  #events to sample for background
  sample.N<-floor(1E4/dt[barcode!=0,length(unique(barcode))])
  #background data
  dt.bkgd<-dt[
    barcode!=0,
    .SD[sample(.N,sample.N,replace = T)],
    keyby=.(barcode),.SDcols = barcode.dims][,!'barcode']
  #metadata for title
  .mdat<-unlist(dt[,lapply(.SD,function(j) unique(as.character(j))),.SDcols = mdat.cols])
  #plot list: split by node
  plist<-sapply(
    split(
      dt[
        barcode==0,
        if(.N>node.threshold){.SD[sample(.N,node.threshold,replace = F)]}else{.SD},
        keyby=.(node_barcode)][
          ,.SD,.SDcols = c(barcode.dims,'node_barcode')],
      by='node_barcode',
      sorted=T
    ),function(dt.node){
      plot.pairs<-lapply(pp,function(.pp){
        ggplot2::ggplot(data=NULL,ggplot2::aes(!!as.name(.pp[1]),!!as.name(.pp[2]))) +
          ggplot2::geom_hex(data=dt.bkgd,fill='gray',bins=bins.N)+
          ggplot2::geom_hex(data=dt.node,bins=bins.N) +
          viridis::scale_fill_viridis(option = "plasma") +
          ggplot2::coord_cartesian(xlim=lims[[.pp[1]]],ylim=lims[[.pp[2]]]) +
          ggplot2::theme_classic() +
          ggplot2::guides(fill = "none")
      })
      patchwork::wrap_plots(plot.pairs,guides = 'collect') +
        patchwork::plot_annotation(title=paste(
          paste("Barcode 0 Node:", dt.node[,unique(node_barcode)]),
          "",
          paste(names(.mdat),.mdat,sep = ": ",collapse = "\n"),
          sep="\n")
        )
    },simplify = F)
  #
  return(plist)
}
#' @title Heatmap of cluster-specific, per-marker median expression values
#' @description
#' A heatmap of median expression values is a critical component of exploring clustering results. This function will return either a static \link[pheatmap]{pheatmap} for displaying/printing or an interactive \link[plotly]{plot_ly} heatmap object.
#'
#' @param fsom An object as returned from \link{som.codes.dt}; fsom$codes.dt is required.
#' @param heatmap.type *Matched argument* Character string used to determine the returned heatmap type: "pheatmap" or "plotly".
#' @param color.function *Matched argument* Character string used to determine which pre-defined heatmap color function to use: 'ryb' or 'greens'. See \link[RColorBrewer]{brewer.pal} **Details**.
#' @param color.override Color function; if defined, will override the default `color.function`.
#' @param color.n Numeric (default = 100); the number of colors to interpolate using the defined `color.function`.
#' @param colnames.split Character/separator (" ","_"); used to split column names to improve 'readability'.
#' @param use.censoring *Specialized case* Logical; if TRUE and a 'censor' column exists in fsom$codes.dt, cluster(s) marked as TRUE will be dropped from the heatmap.
#' @param clusters.to.censor *Specialized case* Numeric vector; the defined cluster(s) will be dropped from the heatmap.
#' @param cols.to.censor *Specialized case* Character vector; the defined column(s) will be dropped from the heatmap.
#' @param ... \link[pheatmap]{pheatmap} arguments. Primarily used to define `scale='row'`.
#'
#' @return A \link[pheatmap]{pheatmap} or \link[plotly]{plot_ly} heatmap, depending on `heatmap.type`. If `heatmap.type='plotly'`, the heatmap object will be row and column ordered to match `pheatmap`s dendrograms.
#' @export
#'
#' @examples
#' fsom<-SOMnambulate:::prepared.examples('fsom')
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=20)
#'
#' ##heatmap; default 'pheatmap'
#' plotSOMheat(fsom)
#'
#' ##heatmap; row-scaled 'pheatmap'
#' plotSOMheat(fsom,scale='row')
#'
#' ##heatmap; row-scaled 'pheatmap'; drop a subset of clusters
#' plotSOMheat(fsom,scale='row',clusters.to.censor=c(1:10))
#'
#' ##heatmap; row-scaled 'pheatmap'; drop a subset of columns
#' plotSOMheat(fsom,scale='row',cols.to.censor=c("CD19_142Nd"))
#'
#' ##heatmap; row-scaled 'pheatmap'; split column names
#' plotSOMheat(fsom,scale='row',colnames.split = "_")
#'
#' ##heatmap; row-scaled 'pheatmap'; split column names; custom colors
#' plotSOMheat(fsom,scale='row',colnames.split = "_",color.override=viridis::plasma)
#'
#' ##heatmap; row-scaled; split column names; plotly object
#' plotSOMheat(fsom,scale='row',colnames.split = "_",heatmap.type='plotly')
#'
#' ##heatmap; row-scaled; split column names; plotly object; modifed plotly object
#' plotly::hide_colorbar(plotSOMheat(fsom,scale='row',colnames.split = "_",heatmap.type='plotly'))
plotSOMheat<-function(fsom,heatmap.type=c("pheatmap","plotly"),color.function=c('ryb','greens'),color.override=NULL,color.n=100,colnames.split=NULL,use.censoring=T,clusters.to.censor=NULL,cols.to.censor=NULL,...){
  ##
  .args<-list(...)
  ##
  if(is.null(fsom$codes.dt)){
    stop("Need fsom as returned from som.codes.dt(...)")
  }
  cluster.col<-grep("cluster",names(fsom$codes.dt),value = T)
  if(use.censoring&'censor' %in% names(fsom$codes.dt)){
    fsom$codes.dt<-fsom$codes.dt[censor==F]
  }
  if(!is.null(clusters.to.censor)){
    fsom$codes.dt<-fsom$codes.dt[!get(cluster.col) %in% clusters.to.censor]
  }
  ##calculate median expression values, per cluster
  dt.medians<-fsom$codes.dt[,
                            j=lapply(.SD,stats::median),
                            .SDcols = colnames(fsom$codes),
                            keyby=cluster.col]
  ##drop columns here
  if(!is.null(cols.to.censor)){
    dt.medians<-dt.medians[,!names(dt.medians) %in% cols.to.censor, with = F]
  }
  ##row names for the heatmap; as.numeric(as.vector(...)) to preserve non-sequential factor ordering
  .rownames<-as.numeric(as.vector(dt.medians[[cluster.col]]))
  ##column names for the heatmap
  .colnames<-get.dt.col.classes(dt.medians)$numeric
  ##split the column names for 'readability'
  if(!is.null(colnames.split)){
    .colnames<-stringr::str_extract(.colnames,paste0("(\\w+)",colnames.split,"(\\w+)"),group = 1)
  }
  ##matrix for heatmap
  mat<-as.matrix(dt.medians[,.SD,.SDcols = is.numeric])
  colnames(mat)<-.colnames
  rownames(mat)<-.rownames
  ##color functions for heatmap
  .colors<-switch(
    match.arg(color.function),
    greens = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name ="Greens")),
    ryb = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,name="RdYlBu")))
  )
  if(!is.null(color.override)){
    .colors<-color.override
  }
  ##interpolate the colors; create a color vector
  .colors<-.colors(color.n)
  ##pheatmap heatmap
  pheat<-pheatmap::pheatmap(
    mat,
    border_color = NA,
    silent=T,
    color = .colors,
    ...
  )
  ##row-scaled matrix if argument (scale = "row") was defined for pheatmap function
  if('scale' %in% names(.args)){
    if(.args$scale=='row'){
      mat<-t(apply(mat,1,function(x){
        (x - mean(x))/stats::sd(x)
      }))
      rownames(mat)<-.rownames
    }
  }
  ##reorder mat based on pheat dendrogram
  mat<-mat[rev(pheat$tree_row$order),pheat$tree_col$order]
  ##plotly heatmap
  plotly.heatmap <- plotly::plot_ly(
    x = colnames(mat),
    y = rownames(mat),
    z = mat,
    type = "heatmap",
    colors = .colors,
    size = 10,
    source = 'axis.selection'
  )
  plotly.heatmap<-plotly::layout(
    p = plotly.heatmap,
    yaxis = list(
      tickfont = list(size = 10),
      type = "category"
    )
  )
  plotly.heatmap <- plotly::layout(
    p = plotly.heatmap,
    xaxis = list(showspikes=T,
                 spikedash="longdash",
                 spikemode="across",
                 spikecolor="purple",
                 spikethickness=1),
    yaxis = list(showspike=T,
                 spikedash="longdash",
                 spikemode="across",
                 spikecolor="orange",
                 spikethickness=1)
  )
  ##return the heatmap
  switch(match.arg(heatmap.type),
         pheatmap = {grid::grid.newpage();return(pheat)},
         plotly = return(plotly.heatmap)
  )
}
#' @title A cluster-faceted bivariate ggplot with a silhouette sampled from all data points
#' @description
#' Bivariate plot pairs are used to explore FlowSOM clustering results by overlaying cluster-specific data onto a silhouette sampled from all available data; this allows for visualizing how cluster-specific expression relates to the data as a whole.
#'
#' @param dt An object as returned from \link{fcs.to.dt} and mapped using \link{map.som.data}.
#' @param ... \link[ggplot2]{aes} arguments (unquoted); essentially x and y. If a `z` aesthetic is provided, the plot will switch to \link[ggplot2]{stat_summary_hex}, using this defined third variable as a 'color-by'.
#' @param bins \link[ggplot2]{geom_hex} argument; numeric vector giving number of bins in both vertical and horizontal directions.
#' @param show.clusters Numeric vector; if defined, the faceted plot will show only those clusters.
#' @param silhouette.n Numeric; the number of events to sample from `dt` to create the silhouette. Set to 2E4 by default.
#'
#' @return a \link[ggplot2]{ggplot} object; faceted by clusters
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples("dt")
#' fsom<-SOMnambulate:::prepared.examples("fsom")
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=20)
#' for(j in colnames(fsom$codes)){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#' map.som.data(fsom,dt)
#' plotSOMclusters(dt,CD4_145Nd,CD8a_162Dy,bins=50)
#' plotSOMclusters(dt,CD4_145Nd,CD8a_162Dy,bins=50,show.clusters=c(4,14))
#' plotSOMclusters(dt,CD4_145Nd,CD8a_162Dy,bins=50,show.clusters=c(4,14)) +
#' ggplot2::theme_classic() + ggplot2::guides(fill='none')
#'
#' plotSOMclusters(dt,CD4_145Nd,CD8a_162Dy,z=CD4_145Nd,bins=50,show.clusters=c(4,14))
plotSOMclusters<-function(dt,...,bins=200,show.clusters=NULL,silhouette.n=2E4){
  ##
  cluster.col = grep("cluster", names(dt), value = T)
  if(length(cluster.col)!=1){
    stop(paste(
      "More than one cluster column found:",
      "subset dt to include only the cluster-column-of-interest; drop the others",
      "Example: dt[,!c('cluster_live')]",
      sep="\n"
    ))
  }
  ##
  sample.sil<-ifelse(dt[,.N]>silhouette.n,T,F)
  ##
  if("z" %in% ...names()){
    ggplot2::ggplot(
      data = if(!is.null(show.clusters)){dt[get(cluster.col) %in% show.clusters]}else{dt},
      ggplot2::aes(...)) +
      ggplot2::stat_summary_hex(bins=bins) +
      viridis::scale_fill_viridis(option="magma") +
      ggplot2::facet_wrap(~get(cluster.col))
  }else{
    ##
    ggplot2::ggplot(data=NULL,ggplot2::aes(...)) +
      ggplot2::geom_hex(
        data=if(dt[,.N]>silhouette.n){
          set.seed(1337)
          dt[sample(.N,silhouette.n),!cluster.col,with=F]
        }else{
          dt[,!cluster.col,with=F]
        },
        bins=bins,fill='gray') +
      ggplot2::geom_hex(
        data=if(!is.null(show.clusters)){dt[get(cluster.col) %in% show.clusters]}else{dt},
        bins=bins
      ) +
      viridis::scale_fill_viridis(
        option="plasma"
        ,limits = c(0,bins=bins)
        ,oob = scales::squish
      ) +
      ggplot2::facet_wrap(~get(cluster.col))
  }
}
#' @title A UMAP of FlowSOM nodes colored by cluster value or expression value
#'
#' @param fsom An object as returned from \link{som.codes.dt}; `fsom$codes.dt` is required.
#' @param dt An object as returned from \link{map.som.data}
#' @param .by Character string; a column in `dt` -- usually a factor. `.by` is used in \link{dt.N.counts} to generate node counts and in \link[ggplot2]{facet_wrap} to generate facets.
#' @param cluster.highlight Numeric; if defined, will highlight specific clusters.
#' @param expression.facets Logical; if `TRUE`, the UMAP will be faceted by variable and colored by (scaled) expression value.
#'
#' @return a \link[ggplot2]{ggplot} object
#' @export
#'
#' @examples
#' dt<-SOMnambulate:::prepared.examples("dt")
#' fsom<-SOMnambulate:::prepared.examples("fsom")
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=20, umap = TRUE)
#' for(j in colnames(fsom$codes)){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#' map.som.data(fsom,dt)
#'
#' plotSOMumap(fsom,dt)
#'
#' #faceted by 'stim.condition'
#' plotSOMumap(fsom,dt,.by='stim.condition')
#'
#' #faceted by 'sample.id'
#' #for datasets with many samples, this is not feasible
#' plotSOMumap(fsom,dt,.by='sample.id')
#'
#' #data.table's 'i' variable to select a single sample
#' plotSOMumap(fsom,dt[sample.id %in% unique(sample.id)[1]],.by='sample.id')
#'
#' p<-plotSOMumap(fsom,dt[sample.id %in% unique(sample.id)[1]],.by='sample.id',cluster.highlight=4)
#' p
#' #modify the ggplot object
#' p + ggplot2::scale_size_area(max_size = 15,guide="none")
#'
#' #colored by marker expression
#' plotSOMumap(fsom,dt,expression.facets=TRUE)
plotSOMumap<-function(fsom,dt,.by=NULL,cluster.highlight=NULL,expression.facets=FALSE){
  if(is.null(fsom$codes.dt)){
    stop("Need fsom as returned from som.codes.dt(...)")
  }
  if(!all(c("umap.1","umap.2") %in% names(fsom$codes.dt))){
    stop("Need fsom umap dimensions as returned from som.codes.dt(...,umap.codes=TRUE)")
  }
  ##
  cluster.col<-grep("cluster",names(fsom$codes.dt),value = T)
  node.col<-sub("cluster","node",cluster.col)
  if(!all(c(cluster.col,node.col) %in% names(dt))){
    stop("Need dt as returned from map.som.data(...)")
  }
  ##
  #dt[,length(unique(get(.by)))]>20#warning conditional for too many facets?
  node.counts<-dt.N.counts(dt[,.N,keyby = c(.by,node.col)])
  #node.sizes<-node.counts[,.(node.size.median=stats::median(N.freq)),keyby = .(node_pbmc)]
  ##
  if(expression.facets){
    #melt, merge and use for plots
    codes.melt<-data.table::melt(
      data=fsom$codes.dt,
      id.vars=grep('node|cluster|umap',names(fsom$codes.dt)))
    codes.melt[,
               value.scaled:= (value-min(value))/(max(value)-min(value)),
               by=variable]
    ggplot2::ggplot(data=merge(codes.melt,node.counts)) +
      ggplot2::aes(umap.1,umap.2,color=value.scaled,size=N.freq) +
      ggplot2::geom_point(alpha=.75) +
      ggplot2::scale_size_area(max_size = 10,guide="none") +
      viridis::scale_color_viridis(option='plasma') +
      ggplot2::theme_dark() +
      ggplot2::theme(panel.grid = ggplot2::element_blank()) +
      ggplot2::facet_wrap(~variable)
  }else{
    ##all clusters
    p<-ggplot2::ggplot(data=merge(fsom$codes.dt,node.counts)) +
      ggplot2::aes(umap.1,umap.2,size=N.freq) +
      ggplot2::theme_void() +
      ggplot2::theme(panel.border = ggplot2::element_rect(color="black", fill=NA)) +
      ggplot2::geom_point(ggplot2::aes(color=!!as.name(cluster.col)),alpha=.5) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3))) +
      ggplot2::facet_wrap(.by)
    ##cluster highlight
    if(!is.null(cluster.highlight)){
      p<- p +
        ggplot2::geom_point(
          data=merge(fsom$codes.dt,node.counts)[get(cluster.col) %in% cluster.highlight],
          color="black",alpha=.5)
    }
    #ggplot2::scale_size_area(max_size = 15,guide="none")
    return(p)
  }
}
#' @title Interactive exploration of high-dimensional clustering (FlowSOM) results
#'
#' @param dt An object as returned from \link{fcs.to.dt} and mapped using \link{map.som.data}.
#' @param fsom An object as returned from \link{som.codes.dt}; fsom$codes.dt is required.
#'
#' @return \link[shiny]{shinyApp}
#' @export
#'
#' @examples
#' if (interactive()) {
#' dt<-SOMnambulate:::prepared.examples("dt")
#' fsom<-SOMnambulate:::prepared.examples("fsom")
#' fsom<-som.codes.dt(fsom,append.name = 'pbmc',k=20,umap.codes=TRUE)
#' for(j in colnames(fsom$codes)){data.table::set(dt,i=NULL,j=j,value=asinh(dt[[j]]/10))}
#' map.som.data(fsom,dt)
#'
#' plotSOMshiny(dt,fsom)
#'}
plotSOMshiny<-function(dt,fsom){
  ##use internal function to get dt column classes
  col.classes<-get.dt.col.classes(dt)
  ##vectors for use as shiny inputs
  choices.numeric<-col.classes$numeric
  choices.factor<-col.classes$factor
  choices.sample<-dt[,sort(unique(sample.id))]
  cluster.col<-grep('cluster',names(fsom$codes.dt),value = T)
  choices.cluster<-levels(fsom$codes.dt[,get(cluster.col)])
  ##limits for numeric columns; all data
  lims <- dt[,lapply(.SD,range),.SDcols=choices.numeric]
  ##limits for UMAP
  lims.umap <- fsom$codes.dt[,lapply(.SD,range),.SDcols=c('umap.1','umap.2')]
  ##cluster counts
  cluster.counts<-dt.N.counts(dt[,.N,keyby=c('sample.id',cluster.col,choices.factor)])
  ##shiny ui inputs
  ui.inputs<-list(
    marker.x=shiny::selectInput(
      inputId = "marker1",
      label = "Marker (x):",
      choices = choices.numeric,
      selected = choices.numeric[1]
    ),
    marker.y=shiny::selectInput(
      inputId = "marker2",
      label = "Marker (y):",
      choices = choices.numeric,
      selected = choices.numeric[2]
    ),
    color.by=shiny::selectInput(
      inputId = "color.by",
      label = "Color:",
      choices = choices.numeric,
      selected = ""
    ),
    color.check=shiny::checkboxInput(
      inputId = "color.check",
      label = "Color by selected variable?",
      value = FALSE
    ),
    sample.id=shiny::selectInput(
      inputId = "sample.id",
      label="Sample:",
      choices=choices.sample,
      selected = choices.sample[1]
    ),
    cluster.val=shiny::selectInput(
      inputId = "cluster.val",
      label="Cluster #:",
      choices=choices.cluster,
      selected = 1
    ),
    factors=shiny::selectInput(
      inputId = "factors",
      label="Factor:",
      choices=choices.factor,
      selected = choices.factor[1]
    )
  )
  ##
  ui <- shinydashboard::dashboardPage(
    header = shinydashboard::dashboardHeader(
      title = "SOMthingShiny",disable = F
    ),
    sidebar = shinydashboard::dashboardSidebar(
      ui.inputs$marker.x,
      ui.inputs$marker.y,
      ui.inputs$color.by,
      ui.inputs$color.check,
      ui.inputs$cluster.val,
      ui.inputs$sample.id,
      ui.inputs$factors
    ),
    body = shinydashboard::dashboardBody(
      shiny::fluidRow(
        ##bivariate plot; silhouette with cluster overlay
        shinydashboard::box(
          collapsible = T,
          title=NULL,
          shiny::plotOutput(".plotSOMclusters"),
          width=5
        ),
        ##UMAP
        shinydashboard::box(
          collapsible = T,
          title=NULL,
          shiny::plotOutput(".plotSOMumap"),
          width=5
        )
      ),
      shiny::fluidRow(
        ##clickable heatmap
        shinydashboard::box(
          collapsible = T,
          title="Axis (X,Y) Selection (click,click)",
          plotly::plotlyOutput("plotly_heat"),
          width = 10
        )
      ),
      shiny::fluidRow(
        ##boxplot and points
        shinydashboard::box(
          collapsible = T,
          title=NULL,
          plotly::plotlyOutput(".plotSOMcounts_plotly")
        )
      )
    )
  )
  ##
  server = function(input, output) {
    ##bivariate plot; silhouette with cluster overlay
    .plotSOMclusters <- shiny::reactive({
      p <- if(input$color.check){
        plotSOMclusters(
          dt[sample.id %in% input$sample.id],
          x=!!as.name(input$marker1),
          y=!!as.name(input$marker2),
          z=!!as.name(input$color.by),
          show.clusters = input$cluster.val
        )
      }else{
        plotSOMclusters(
          dt[sample.id %in% input$sample.id],
          x=!!as.name(input$marker1),
          y=!!as.name(input$marker2),
          show.clusters = input$cluster.val
        )
      }
      p <- p +
        ggplot2::theme_classic() +
        ggplot2::guides(fill='none') +
        if(all(c(input$marker1,input$marker2) %in% names(lims))){
          ggplot2::coord_cartesian(
            xlim=lims[,get(input$marker1)],
            ylim=lims[,get(input$marker2)])
        }else if(input$marker1 %in% names(lims)){
          ggplot2::coord_cartesian(
            xlim=lims[,get(input$marker1)])
        }else if(input$marker2 %in% names(lims)){
          ggplot2::coord_cartesian(
            ylim=lims[,get(input$marker2)])
        }
      return(p)
    })
    ##bivariate plot; output
    output$.plotSOMclusters <- shiny::renderPlot({
      .plotSOMclusters()
    })
    ##umap plot
    .plotSOMumap <- shiny::reactive({
      plotSOMumap(
        fsom,
        dt[sample.id %in% input$sample.id],
        .by='sample.id',
        cluster.highlight = input$cluster.val
      ) +
        ggplot2::scale_size_area(guide="none") +
        ggplot2::coord_cartesian(
          xlim=lims.umap[['umap.1']],
          ylim=lims.umap[['umap.2']]
        )
    })
    ##umap plot; output
    output$.plotSOMumap <- shiny::renderPlot({
      .plotSOMumap()
    })
    ##clickable heatmap; data.frame to 'hold' click values
    clicks <- shiny::reactiveValues(dat = data.frame(marker1 = NA, marker2 = NA))
    ##clickable heatmap; reactive
    click <- shiny::reactive({
      plotly::event_data("plotly_click", priority = 'event', source = 'axis.selection')
    })
    ##clickable heatmap; update data.frame with reactive 'clicks'
    shiny::observeEvent(eventExpr = click(),{
      if(is.na(clicks$dat$marker1)&is.na(clicks$dat$marker2)){
        clicks$dat$marker1 <- click()$x
      }else if(is.na(clicks$dat$marker2)){
        clicks$dat$marker2 <- click()$x
        shiny::updateSelectInput(inputId = 'marker1',
                                 selected = clicks$dat$marker1)
        shiny::updateSelectInput(inputId = 'marker2',
                                 selected = clicks$dat$marker2)
        shiny::updateSelectInput(inputId = "cluster.val",
                                 selected = stringr::str_extract(click()$y,
                                                                 "[0-9]+"))
        clicks$dat$marker1 <- NA
        clicks$dat$marker2 <- NA
      }
    })
    ##clickable heatmap; output
    output$plotly_heat <- plotly::renderPlotly(
      plotly::hide_colorbar(
        plotSOMheat(
          fsom,
          heatmap.type = 'plotly',
          scale="row"
        )
      )
    )
    ##boxplot with points; subset by cluster
    .plotSOMcounts_gg <- shiny::reactive({
      ggplot2::ggplot(cluster.counts[get(cluster.col)==input$cluster.val],ggplot2::aes(!!as.name(input$factors),N.freq,label.1=sample.id)) +
        ggplot2::geom_boxplot(,outlier.shape = NA) +
        ggplot2::geom_jitter(width=0.25,alpha=0.2) +
        ggplot2::facet_wrap(cluster.col)
    })
    ##boxplot with points; subset by cluster; output -- ggplot
    output$.plotSOMcounts_gg <- shiny::renderPlot({
      .plotSOMcounts_gg()
    })
    ##boxplot with points; subset by cluster; output -- plotly
    output$.plotSOMcounts_plotly <- plotly::renderPlotly({
      p<-plotly::ggplotly(.plotSOMcounts_gg(),source="boxplot.click")
      p$x$data[[1]]$marker$opacity<-0
      return(p)
    })
    ##reactive boxplot clicks
    click.boxplot <- shiny::reactive({
      plotly::event_data("plotly_click", priority = 'event', source = 'boxplot.click')
    })
    shiny::observeEvent(eventExpr = click.boxplot(),{
      #print(click.boxplot())
      .number<-click.boxplot()$pointNumber+1
      .val<-click.boxplot()$y
      .sample.id<-cluster.counts[get(cluster.col)==input$cluster.val,sample.id[.number]]
      shiny::updateSelectInput(
        inputId = "sample.id",
        selected = .sample.id
      )
    })
  }
  ##
  shiny::shinyApp(ui, server)
}
