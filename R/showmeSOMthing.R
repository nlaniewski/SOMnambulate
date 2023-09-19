#' @title Interactive exploration of high-dimensional clustering (FlowSOM) results
#' @description
#' This interactive \code{shiny}\\code{shinydashboard} allows a visual exploration of FlowSOM clustering results. \code{plotly} clickable cluster heatmap dimensions and factored datapoints are the primary features.
#'
#'
#' @param dat \code{data.table} containing a factored \code{FlowSOM} cluster column
#' @param code.medians pre-calculated cluster medians; a listed result of \code{fsom.codes.to.clusters}
#' @param marker.pair character vector of length 2; if defined, will pre-select the initially displayed dimensions for the bivariate plots.  The default \code{NULL} will pre-select the first two columns of dat for display.
#' @param asinh.view Logical. By default, \code{FALSE}; if \code{TRUE}, will display an interactive plot to assess the effect of \code{ashinh} transformations on the primary data.
#' @param use.cluster.counts Logical. By default, \code{TRUE}; will key \code{dat} by the unique identifier 'sample.id' and cluster column to generate per-sample cluster counts; cluster counts are used as input for factored boxplots.
#'
#' @return \code{shinyapp}
#' @export
#'
showmeSOMething <- function(dat,code.medians=NULL,marker.pair=NULL,asinh.view=F,use.cluster.counts=T){
  if(!data.table::is.data.table(dat)) stop("Need a data.table...")
  ##
  dat.names <- names(dat)
  cluster.col<-cluster.col.check(dat.names)
  ##
  if(!all(c('sample.id',cluster.col) %in% names(dat))) stop("Need 'sample.id' and 'pbmc_cluster' columns in the data.table")
  ##
  if(use.cluster.counts){
    data.table::setkeyv(dat,c('sample.id',cluster.col))
    cluster.counts<-dat[,unique(.SD),.SDcols = !is.numeric][dat[,.N,keyby = c('sample.id',cluster.col)]]
    cluster.counts[,'percent.of' := N/sum(N)*100,by='sample.id']
    ##
    clusters.seq<-as.numeric(unique(cluster.counts[[cluster.col]]))
    ##
    samples <- unique(cluster.counts[['sample.id']])
    ##
    factor.choices <- grep(cluster.col,names(which(sapply(cluster.counts,is.factor))),invert = T,value = T)
    value.choices <- names(which(sapply(cluster.counts,is.numeric)))
  }else{
    samples <- unique(dat[['sample.id']])
  }
  ##
  if(!is.null(marker.pair)){
    m1 <- marker.pair[1];m2 <- marker.pair[2]
  }else{
    if(all(c('FSC_A','SSC_A') %in% dat.names)){
      m1 <- 'FSC_A'
      m2 <- 'SSC_A'
    }else{
      m1 <- dat.names[1]
      m2 <- dat.names[2]
    }
  }
  ##
  lims <- dat[,lapply(.SD,range),.SDcols=is.numeric]
  ##
  axis.click.select<-shiny::fluidRow(
    shinydashboard::box(
      title="Axis (X,Y) Selection",
      plotly::plotlyOutput("plotly_heat",height="200px"),
      width=10
    )
  )
  ##
  if(asinh.view==T){
    slider.vals <- stats::setNames(nm=c('min','max','value','step'),
                                   c(1,3000,1000,50)
    )
    # any(dat==0)&dat[,.N]*ncol(dat)/length(which(dat==0))>1
    if(sum(dat[,lapply(.SD,function(x)length(which(x==0)))])/dat[,.N]*ncol(dat)*100>1){#preponderance of zeroes in mass cyto. data
      slider.vals[1:4] <- c(1,10,2,1)
    }
  }
  ##shinydashboard items;store as variables
  par.menu <- shinydashboard::menuItem(
    "Parameter Selection:",
    tabName = "markers",
    shiny::selectInput(inputId = "marker1",
                       label = "Marker (x):",
                       choices = dat.names,
                       selected = m1),
    shiny::selectInput(inputId = "marker2",
                       label = "Marker (y):",
                       choices = dat.names,
                       selected = m2),
    shiny::numericInput(inputId = "rowsamp",
                        label = "# of 'Events' to display:",
                        value = 1E4,
                        min = 1E4,
                        max = 1E4,
                        step = 1E4),
    shiny::selectInput(inputId = "sample.id",
                       label="Sample",
                       choices=samples,
                       selected = samples[1]),
    if(!is.null(clusters.seq)){
      shiny::selectInput(inputId = "cluster.val",
                         label="Cluster #",
                         choices=clusters.seq,
                         selected = NULL)
    },
    shiny::radioButtons(inputId = "axis.select",
                        label = "Axis Select Type:",
                        choices=if(!is.null(clusters.seq)) c("Markers","Clusters") else "Markers",
                        selected = if(!is.null(clusters.seq)) "Clusters" else "Markers",
                        inline = T)
  )
  factor.menu <- shinydashboard::renderMenu({
    shinydashboard::menuItem(
      "Factor Selection:",
      tabName = "factors",
      shiny::selectInput(inputId = "factor.name",
                         label = "Factor (x):",
                         choices = unname(factor.choices),
                         selected = 'visit'),
      shiny::radioButtons(inputId = "value.y",
                          label = "Value (y):",
                          choices= value.choices,
                          selected = "prop",
                          inline = T)
    )
  })
  asinh.menu <- shinydashboard::renderMenu({
    shinydashboard::menuItem(
      "Asinh Transform:",
      tabName = "asinh",
      shiny::radioButtons(
        inputId = 'asinh.applied',
        label="Apply asinh?",
        choices = c('Yes',"No"),
        selected = 'No',
        inline = T
      ),
      shiny::sliderInput(
        inputId = 'cofactor.xaxis',
        label = "Cofactor: X-axis",
        min = slider.vals[['min']],
        max = slider.vals[['max']],
        value = slider.vals[['value']],
        step = slider.vals[['step']]
      ),
      shiny::sliderInput(
        inputId = 'cofactor.yaxis',
        label = "Cofactor: Y-axis",
        min = slider.vals[['min']],
        max = slider.vals[['max']],
        value = slider.vals[['value']],
        step = slider.vals[['step']]
      )
    )
  })
  cyto.plots<-shiny::fluidRow(
    shinydashboard::box(
      collapsible = T,
      title=NULL,
      shiny::plotOutput("ggbivariate_plot1"),
      width=5
    ),
    shinydashboard::box(
      collapsible = T,
      collapsed = is.null(clusters.seq),
      title=NULL,
      shiny::plotOutput("ggbivariate_plot2"),
      width=5
    )
  )
  axis.click.select<-shiny::fluidRow(
    shinydashboard::box(
      collapsible = T,
      title="Axis (X,Y) Selection",
      plotly::plotlyOutput("plotly_heat"),
      width = 10
    )
  )
  factor.plot<-shiny::fluidRow(
    shinydashboard::box(
      collapsible = T,
      collapsed = is.null(clusters.seq),
      title=NULL,
      plotly::plotlyOutput("factor_plot1"),
      width=10
    )
  )
  ##
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = paste("showmeSOMthing"),
                                    disable = F),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        par.menu,
        shinydashboard::menuItemOutput('factor.menu'),
        shinydashboard::menuItemOutput('asinh.menu')
      )
    ),
    shinydashboard::dashboardBody(
      cyto.plots,
      axis.click.select,
      factor.plot
    )
  )

  server <- function(input, output) {
    ##
    if(asinh.view==TRUE){
      output$asinh.menu <- asinh.menu
    }
    ##
    if(!is.null(clusters.seq)){
      output$factor.menu <- factor.menu
    }
    ##
    shiny::observeEvent(input$sample.id,{
      shiny::updateNumericInput(inputId = "rowsamp",
                                value = ifelse(dat[get('sample.id')==input$sample.id ,.N]<1E5,dat[get('sample.id')==input$sample.id ,.N],1E5),
                                min = 1E5,
                                max = dat[get('sample.id')==input$sample.id ,.N],
                                step = 1E5)
    })
    ##
    ggbivariate_plot1 <- shiny::reactive({
      if(asinh.view==FALSE){
        p.tmp <- gg.func.bivariate(dat = if(dat[get('sample.id')==input$sample.id,.N]>input$rowsamp){
          dat[get('sample.id')==input$sample.id][sample(.N,input$rowsamp)]
        }else{
          dat[get('sample.id')==input$sample.id]
        },
        x = !!ggplot2::sym(input$marker1),
        y = !!ggplot2::sym(input$marker2)
        )
      }else if(asinh.view==TRUE){
        shiny::req(input$asinh.applied)
        if(input$asinh.applied=='Yes'){
          p.tmp <- gg.func.bivariate(dat = if(dat[get('sample.id')==input$sample.id,.N]>input$rowsamp){
            dat[get('sample.id')==input$sample.id][sample(.N,input$rowsamp)]
          }else{
            dat[get('sample.id')==input$sample.id]
          },
          x = asinh(!!ggplot2::sym(input$marker1)/input$cofactor.xaxis),
          y = asinh(!!ggplot2::sym(input$marker2)/input$cofactor.yaxis)
          )
        }else if(input$asinh.applied=='No'){
          p.tmp <- gg.func.bivariate(dat = if(dat[get('sample.id')==input$sample.id,.N]>input$rowsamp){
            dat[get('sample.id')==input$sample.id][sample(.N,input$rowsamp)]
          }else{
            dat[get('sample.id')==input$sample.id]
          },
          x = !!ggplot2::sym(input$marker1),
          y = !!ggplot2::sym(input$marker2)
          )
        }
      }
      p.tmp <- p.tmp +
        ggplot2::labs(title = "All Events",
                      subtitle = paste(input$rowsamp, "of", dat[get('sample.id')==input$sample.id,.N], "displayed")) +
        ggplot2::xlab(input$marker1) +
        ggplot2::ylab(input$marker2)
      if(asinh.view==TRUE){
        return(p.tmp)
      }else{
        if(all(c(input$marker1,input$marker2) %in% names(lims))){
          p.tmp <- p.tmp +
            ggplot2::coord_cartesian(xlim=lims[,get(input$marker1)],
                                     ylim=lims[,get(input$marker2)])
        }else if(input$marker1 %in% names(lims)){
          p.tmp <- p.tmp +
            ggplot2::coord_cartesian(xlim=lims[,get(input$marker1)])
        }else if(input$marker2 %in% names(lims)){
          p.tmp <- p.tmp +
            ggplot2::coord_cartesian(ylim=lims[,get(input$marker2)])
        }
        return(p.tmp)
      }
    })
    ##
    if(!is.null(clusters.seq)){
      ggbivariate_plot2 <- shiny::reactive({
        p.tmp <- gg.func.bivariate(dat[get('sample.id')==input$sample.id][get(cluster.col)==input$cluster.val],
                                   x = !!ggplot2::sym(input$marker1),
                                   y = !!ggplot2::sym(input$marker2))
        p.tmp <- p.tmp +
          ggplot2::labs(title = paste("Cluster #",input$cluster.val),
                        subtitle = paste(dat[get('sample.id')==input$sample.id][get(cluster.col)==input$cluster.val,.N], "displayed")) +
          ggplot2::xlab(input$marker1) +
          ggplot2::ylab(input$marker2)
        if(all(c(input$marker1,input$marker2) %in% names(lims))){
          p.tmp <- p.tmp +
            ggplot2::coord_cartesian(xlim=lims[,get(input$marker1)],
                                     ylim=lims[,get(input$marker2)])
        }else if(input$marker1 %in% names(lims)){
          p.tmp <- p.tmp +
            ggplot2::coord_cartesian(xlim=lims[,get(input$marker1)])
        }else if(input$marker2 %in% names(lims)){
          p.tmp <- p.tmp +
            ggplot2::coord_cartesian(ylim=lims[,get(input$marker2)])
        }
        return(p.tmp)
      })
    }
    ##
    #if(!is.null(clusters.seq)){
    factor_plot1 <- shiny::reactive({
      shiny::req(input$factor.name,input$value.y)
      p<-plotly::ggplotly(
        gg.func.boxplot.points(cluster.counts[get(cluster.col)==input$cluster.val],
                               x=!!ggplot2::sym(input$factor.name),
                               y=!!ggplot2::sym(input$value.y)
        ) + ggplot2::labs(y=if(input$value.y=='percent.of'){
          paste("% of",paste0(toupper(stringr::str_extract(cluster.col,"^[^_]+")),"+"),"Clusters")
        }else if(input$value.y=='N'){
          paste("Number of cells per",paste0(toupper(stringr::str_extract(cluster.col,"^[^_]+")),"+"),"Clusters")
        }),
        tooltip = 'text',
        source = 'sample.selection'
      )
      for(i in factor.choices){
        p$x$data[[2]]$text <- paste0(p$x$data[[2]]$text,'<br />',
                                     paste(paste0(toupper(i),":"),cluster.counts[get(cluster.col)==input$cluster.val][[i]]))
      }
      return(p)
    })
    #}
    ##
    output$ggbivariate_plot1 <- shiny::renderPlot({
      ggbivariate_plot1()
    })
    if(!is.null(clusters.seq)){
      output$ggbivariate_plot2 <- shiny::renderPlot({
        ggbivariate_plot2()
      })
    }else{
      output$ggbivariate_plot2 <- shiny::renderPlot({
        NULL
      })
    }
    ##
    shiny::observeEvent(input$axis.select,{
      if(input$axis.select=="Markers"){
        output$plotly_heat <- plotly::renderPlotly(axis.selection.plotly.heatmap(dat.names))
      }else if(input$axis.select=="Clusters"){
        output$plotly_heat <- plotly::renderPlotly(cluster.axis.selection.plotly.heatmap(code.medians))
      }
    })
    ##
    clicks <- shiny::reactiveValues(dat = data.frame(marker1 = NA, marker2 = NA))
    click <- shiny::reactive({
      plotly::event_data("plotly_click", priority = 'event', source = 'axis.selection')
    })
    ##
    #if(!is.null(clusters.seq)){
    sample_click <- shiny::reactive({
      shiny::req(input$factor.name,input$value.y)
      plotly::event_data("plotly_click", priority = 'event', source = 'sample.selection')
    })
    #}
    ##
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
    ##
    #if(!is.null(clusters.seq)){
    shiny::observeEvent(eventExpr = sample_click(),{
      pn <- sample_click()$pointNumber+1
      shiny::updateSelectInput(inputId = "sample.id",
                               selected = cluster.counts[,unique(get('sample.id'))][pn]
      )
    })
    #}
    ##
    if(!is.null(clusters.seq)){
      output$factor_plot1 <- plotly::renderPlotly(factor_plot1())
    }else{
      output$factor_plot1 <- plotly::renderPlotly(NULL)
    }
    ##
  }
  ##
  shiny::shinyApp(ui, server)
}

gg.func.bivariate <- function(dat,...,bins=100,fill.limits=c(0,50)){
  ggplot2::ggplot(dat,ggplot2::aes(...)) +
    ggplot2::geom_hex(bins = bins) +
    viridis::scale_fill_viridis(option = "plasma", limits = fill.limits, oob = scales::squish)
}

gg.func.boxplot.points <- function(dat,...){
  cluster.col<-cluster.col.check(names(dat))
  ggplot2::ggplot(dat,ggplot2::aes(...)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(size=1.5) +
    ggplot2::facet_wrap(cluster.col)
}

axis.selection.plotly.heatmap<-function(column.names){
  axis.selection<-matrix(1:length(column.names),
                         nrow=1,
                         ncol=length(column.names),
                         dimnames = list(c('click.x_click.y'),column.names)
  )

  plotly.heatmap <- plotly::plot_ly(x=colnames(axis.selection),
                                    y=rownames(axis.selection),
                                    z=axis.selection,
                                    zauto = F,
                                    type = "heatmap",
                                    xgap = 2,
                                    text=colnames(axis.selection),
                                    source = 'axis.selection'
  )
  plotly.heatmap <- plotly::layout(plotly.heatmap,
                                   yaxis = list(tickfont = list(size = 10),
                                                type = "category")
  )
  plotly.heatmap <- plotly::layout(plotly.heatmap,
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
  plotly.heatmap<-plotly::hide_colorbar(plotly.heatmap)
}

cluster.axis.selection.plotly.heatmap<- function(code.medians,row.scale=T,break.vec=NULL,color.type=c("pheatmap.default","greens")){
  if(data.table::is.data.table(code.medians)){
    m<-as.matrix(code.medians[,.SD,.SDcols = is.numeric]);m.cols<-colnames(m)
  }else{
    stop("Need a data.table; median values of FSOM::SOM$codes, grouped by 'cluster'")
  }
  ##row-scale
  if(row.scale){
    m<-t(apply(m,1,scale));colnames(m)<-m.cols
  }
  #prepare cluster ordering for plotly; uses dendgrogram sorting
  row.order <- rev(unlist(stats::as.dendrogram(stats::hclust(stats::dist(m)))))
  col.order <- unlist(stats::as.dendrogram(stats::hclust(stats::dist(t(m)))))
  m<-m[row.order,col.order]
  rownames(m)<-row.order
  ##
  if(is.null(break.vec)){
    r<-range(apply(m,2,range))
    break.vec<-seq(floor(min(r)),ceiling(max(r)),by=0.05)
  }
  ##
  #plotly heatmap
  plotly.heatmap <- plotly::plot_ly(x=colnames(m),
                                    y=rownames(m),
                                    z=m,
                                    zauto = F,
                                    zmin = min(break.vec),
                                    zmax = max(break.vec),
                                    type = "heatmap",
                                    colors = switch(match.arg(color.type),
                                                    pheatmap.default = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,name="RdYlBu")))(length(break.vec)),
                                                    greens = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name ="Greens"))(length(break.vec))
                                    ),
                                    source = 'axis.selection'
  )
  plotly.heatmap <- plotly::layout(plotly.heatmap,
                                   yaxis = list(tickfont = list(size = 10),
                                                type = "category")
  )
  plotly.heatmap <- plotly::layout(plotly.heatmap,
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
  plotly.heatmap<-plotly::hide_colorbar(plotly.heatmap)
  return(plotly.heatmap)
}

gg.func.bivariate.cluster.overlay <- function(dat,...,bins=100,fill.limits=c(0,50),cluster.number=NULL,
                                              overlay.total=1E5,per.sample.cluster.total=5E4,use.labs=F){
  cluster.col<-cluster.col.check(names(dat))
  p<-ggplot2::ggplot(dat[sample(.N,overlay.total)],ggplot2::aes(...)) +
    ggplot2::geom_hex(fill = "gray", bins = bins) +
    ggplot2::geom_hex(data=subsample_dt_sample.id(dat[get(cluster.col)==cluster.number],per.sample.cluster.total),bins=bins) +
    viridis::scale_fill_viridis(option = "plasma", limits = fill.limits, oob = scales::squish) +
    ggplot2::theme_classic() +
    ggplot2::guides(fill='none')
  if(use.labs){
    p<-p+ggplot2::labs(title=paste("Cluster #:",cluster.number),
                       subtitle = paste("Cell Type:",unique(dat[get(cluster.col)==cluster.number,cell.type])),
                       caption = paste(paste("# of events (cluster):",dat[get(cluster.col)==cluster.number,.N]),
                                       paste("# of events (overlay):",overlay.total),
                                       paste('# of events (total)',dat[,.N]),
                                       sep="\n")
    )
  }
  return(p)
}

