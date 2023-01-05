cytoplot <- function(dat,marker.pair=NULL,asinh.view=F){
  if(!data.table::is.data.table(dat)) stop("Need a data.table returned from 'readFCS_dt'...")
  ##
  c.names <- names(dat)
  ##
  if(!is.null(marker.pair)){
    m1 <- marker.pair[1];m2 <- marker.pair[2]
  }else{
    if(all(c('FSC-A','SSC-A') %in% c.names)){
      m1 <- 'FSC-A'
      m2 <- 'SSC-A'
    }else{
      m1 <- c.names[1]
      m2 <- c.names[2]
    }
  }
  ##
  samples<-dat[,unique(sample)]
  ##
  lims <- dat[,lapply(.SD,function(x) c(min(x),max(x))),.SDcols=is.numeric]
  ##
  if('cluster' %in% names(dat)){
    clusters<-dat[,sort(unique(cluster))]
    dat.N.cluster<-cluster_counts_long(dat)
  }else{
    clusters<-NULL
    message("'clusters' is NULL")
  }
  ##
  # if('cell.type' %in% names(dat)){
  #   cell.types<-dat[,unique(na.omit(cell.type))]
  # }
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
    if(any(dat==0)&dat[,.N]*ncol(dat)/length(which(dat==0))>1){#preponderance of zeroes in mass cyto. data
      slider.vals[1:4] <- c(1,10,2,1)
    }
  }
  ##shinydashboard items;store as variables
  par.menu <- shinydashboard::menuItem(
    "Parameter Selection:",
    tabName = "markers",
    shiny::selectInput(inputId = "marker1",
                       label = "Marker (x):",
                       choices = c.names,
                       selected = m1),
    shiny::selectInput(inputId = "marker2",
                       label = "Marker (y):",
                       choices = c.names,
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
    if(!is.null(clusters)){
      shiny::selectInput(inputId = "cluster.val",
                         label="Cluster #",
                         choices=clusters,
                         selected = NULL)
    },
    if(!is.null(clusters)){
      shiny::radioButtons(inputId = "axis.select",
                          label = "Axis Select Type:",
                          choices= c("Markers","Clusters"),
                          selected = "Markers",
                          inline = T)
    }else{
      shiny::radioButtons(inputId = "axis.select",
                          label = "Axis Select Type:",
                          choices= c("Markers"),
                          selected = "Markers",
                          inline = T)
    }
  )
  factor.menu <- shinydashboard::menuItem(
    "Factor Selection:",
    tabName = "factors",
    shiny::selectInput(inputId = "factor.name",
                       label = "Factor (x):",
                       choices = c('subject','visit','condition','batch'),
                       selected = 'visit'),
    shiny::radioButtons(inputId = "value.y",
                        label = "Value (y):",
                        choices= c("prop","per1million"),
                        selected = "prop",
                        inline = T)
  )
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
      title=NULL,
      plotly::plotlyOutput("factor_plot1"),
      width=10
    )
  )
  ##
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = paste("Cyto Plot"),
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
    if(!is.null(clusters)){
      output$factor.menu <- factor.menu
    }
    ##
    shiny::observeEvent(input$sample.id,{
      shiny::updateNumericInput(inputId = "rowsamp",
                                value = ifelse(dat[sample==input$sample.id ,.N]<1E5,dat[sample==input$sample.id ,.N],1E5),
                                min = 1E5,
                                max = dat[sample==input$sample.id ,.N],
                                step = 1E5)
    })
    ##
    ggbivariate_plot1 <- shiny::reactive({
      if(asinh.view==FALSE){
        p.tmp <- gg.func.bivariate(dat[sample==input$sample.id][sample(.N,input$rowsamp)],
                                   x = !!ggplot2::sym(input$marker1),
                                   y = !!ggplot2::sym(input$marker2))
      }else if(asinh.view==TRUE){
        shiny::req(input$asinh.applied)
        if(input$asinh.applied=='Yes'){
          p.tmp <- gg.func.bivariate(dat[sample==input$sample.id][sample(.N,input$rowsamp)],
                                     x = asinh(!!ggplot2::sym(input$marker1)/input$cofactor.xaxis),
                                     y = asinh(!!ggplot2::sym(input$marker2)/input$cofactor.yaxis)
          )
        }else if(input$asinh.applied=='No'){
          p.tmp <- gg.func.bivariate(dat[sample==input$sample.id][sample(.N,input$rowsamp)],
                                     x = !!ggplot2::sym(input$marker1),
                                     y = !!ggplot2::sym(input$marker2)
          )
        }
      }
      p.tmp <- p.tmp +
        ggplot2::labs(title = "All Events",
                      subtitle = paste(input$rowsamp, "of", dat[sample==input$sample.id,.N], "displayed")) +
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
    if(!is.null(clusters)){
      ggbivariate_plot2 <- shiny::reactive({
        p.tmp <- gg.func.bivariate(dat[sample==input$sample.id][cluster==input$cluster.val],
                                   x = !!ggplot2::sym(input$marker1),
                                   y = !!ggplot2::sym(input$marker2))
        p.tmp <- p.tmp +
          ggplot2::labs(title = paste("Cluster #",input$cluster.val),
                        subtitle = paste(dat[sample==input$sample.id][cluster==input$cluster.val,.N], "displayed")) +
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
    if(!is.null(clusters)){
      factor_plot1 <- shiny::reactive({
        plotly::ggplotly(
          gg.func.boxplot.points(dat.N.cluster[cluster==input$cluster.val],
                                 x=!!ggplot2::sym(input$factor.name),
                                 y=!!ggplot2::sym(input$value.y)
          ),
          tooltip = 'text',
          source = 'sample.selection')
      })
    }
    ##
    output$ggbivariate_plot1 <- shiny::renderPlot({
      ggbivariate_plot1()
    })
    if(!is.null(clusters)){
      output$ggbivariate_plot2 <- shiny::renderPlot({
        ggbivariate_plot2()
      })
    }else{
      output$ggbivariate_plot2 <- shiny::renderPlot({
        NULL
      })
    }
    ##
    shiny::observeEvent(input$sample.id,{
      shiny::observeEvent(input$axis.select,{
        if(input$axis.select=="Markers"){
          output$plotly_heat <- plotly::renderPlotly(axis.selection.plotly.heatmap(c.names))
        }else if(input$axis.select=="Clusters"){
          cm=generate_cluster_medians(dat[sample==input$sample.id])[,!'sample']
          output$plotly_heat <- plotly::renderPlotly(cluster.axis.selection.plotly.heatmap(cm))
        }
      })
    })
    ##
    clicks <- shiny::reactiveValues(dat = data.frame(marker1 = NA, marker2 = NA))
    click <- shiny::reactive({
      plotly::event_data("plotly_click", priority = 'event', source = 'axis.selection')
    })
    ##
    if(!is.null(clusters)){
    sample_click <- shiny::reactive({
      plotly::event_data("plotly_click", priority = 'event', source = 'sample.selection')
    })
    }
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
    if(!is.null(clusters)){
    shiny::observeEvent(eventExpr = sample_click(),{
      pn <- sample_click()$pointNumber+1
      shiny::updateSelectInput(inputId = "sample.id",
                               selected = dat.N.cluster[,unique(sample)][pn]
      )
    })
    }
    ##
    if(!is.null(clusters)){
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
  ggplot2::ggplot(dat,ggplot2::aes(...,
                                   text=paste("Subject:",get('subject'),
                                              "<br>Visit:",get('visit'),
                                              "<br>Condition:",get('condition'),
                                              "<br>Batch:",get('batch')
                                   ))) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(size=1.5) +
      ggplot2::facet_wrap(~cluster)
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

cluster.axis.selection.plotly.heatmap<- function(cluster.medians,use.sorting=T,break.vec = seq(0, 1, by = 0.05)){
  if(!is.matrix(cluster.medians)){
    cluster.medians <- as.matrix(cluster.medians)
  }
  #prepare cluster ordering for plotly; uses dendgrogram sorting
  if(use.sorting){
    c.order <- rev(unlist(stats::as.dendrogram(stats::hclust(stats::dist(cluster.medians)))))
    m.order <- unlist(stats::as.dendrogram(stats::hclust(stats::dist(t(cluster.medians)))))
    cluster.medians<-cluster.medians[c.order,m.order]
    rownames(cluster.medians)<-c.order
  }else{
    rownames(cluster.medians)<-seq(nrow(cluster.medians))
  }
  ##
  color.breaks<-grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name ="Greens"))(length(break.vec))
  #plotly heatmap
  plotly.heatmap <- plotly::plot_ly(x=colnames(cluster.medians),
                                    y=rownames(cluster.medians),
                                    z=cluster.medians,
                                    zauto = F,
                                    zmin = min(break.vec),
                                    zmax = max(break.vec),
                                    type = "heatmap",
                                    colors = color.breaks,
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

generate_cluster_medians<-function(dat,use.scale.func=T,by.factor='sample'){
  if(!data.table::is.data.table(dat)&'cluster' %in% names(dat)){
    stop("Need a data.table with a 'cluster' column")
  }
  ##
  cols.for.cluster.medians <- names(which(sapply(dat,is.numeric)[!sapply(dat,is.integer)]))
  cols.for.cluster.medians<-cols.for.cluster.medians[!cols.for.cluster.medians %in% "Time"]
  if(is.null(by.factor)){
    cluster.medians<-dat[,lapply(.SD,stats::median),keyby='cluster',.SDcols=cols.for.cluster.medians][,!'cluster']
  }else{
    cluster.medians<-dat[,lapply(.SD,stats::median),keyby=c('cluster',by.factor),.SDcols=cols.for.cluster.medians][,!'cluster']
  }
  if(use.scale.func){
    if(is.null(by.factor)){
      cluster.medians[,(cols.for.cluster.medians):=lapply(.SD,function(x)(x-mean(x))/stats::sd(x))]
    }else{
      cluster.medians[,(cols.for.cluster.medians):=lapply(.SD,function(x)(x-mean(x))/stats::sd(x)),by=by.factor]
    }
  }
  return(cluster.medians)
}

gg.func.bivariate.cluster.overlay <- function(dat,...,bins=100,fill.limits=c(0,50),cluster.number=NULL,
                                              overlay.total=1E5,use.labs=T){
  p<-ggplot2::ggplot(dat[sample(.N,overlay.total)],ggplot2::aes(...)) +
    ggplot2::geom_hex(fill = "gray", bins = bins) +
    ggplot2::geom_hex(data=dat[cluster==cluster.number],bins=bins) +
    viridis::scale_fill_viridis(option = "plasma", limits = fill.limits, oob = scales::squish) +
    ggplot2::theme_classic() +
    ggplot2::guides(fill='none')
  if(use.labs){
    p<-p+ggplot2::labs(title=paste("Cluster #:",cluster.number),
                       subtitle = paste("Cell Type:",unique(dat[cluster==cluster.number,cell.type])),
                       caption = paste(paste("# of events (cluster):",dat[cluster==cluster.number,.N]),
                                       paste("# of events (overlay):",overlay.total),
                                       paste('# of events (total)',dat[,.N]),
                                       sep="\n")
    )
  }
  return(p)
}
