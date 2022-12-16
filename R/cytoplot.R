cytoplot <- function(dat=NULL,fcs.file.path=NULL,marker.pair=NULL,asinh.view=F){
  if(!is.null(fcs.file.path)){
    dat<-readFCS_dt(fcs.file.path)
  }
  ##
  c.names <- names(dat)
  ##
  if(!is.null(marker.pair)){
    m1 <- marker.pair[1];m2 <- marker.pair[2]
  }else{
    if(all(c('FSC-A','SSC-A') %in% c.names)){
      m1 <- 'FSC-A'
      m2 <- 'SSC-A'
    }
  }
  ##
  if('cluster' %in% names(dat)){
    cm<-generate_cluster_medians(dat)
    p <- cluster.axis.selection.plotly.heatmap(cm)
    axis.click.select<-shiny::fluidRow(
      shinydashboard::box(
        title="Cluster (X,Y) Selection",
        plotly::plotlyOutput("plotly_heat")
      )
    )
  }else{
    p <- axis.selection.plotly.heatmap(c.names)
    axis.click.select<-shiny::fluidRow(
      shinydashboard::box(
        title="Axis (X,Y) Selection",
        plotly::plotlyOutput("plotly_heat",height="150px"),
        width=8
      )
    )
  }
  ##
  total.rows <- nrow(dat)
  if(total.rows<1E5){
    rows.initial<-total.rows
  }else{
    rows.initial<-1E5
  }
  ##shinydashboard items;store as variables
  marker.menu <- shinydashboard::menuItem(
    "Marker Selection:",
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
                        label = "# of 'All Events' to display:",
                        value = rows.initial,
                        min = 1E5,
                        max = total.rows,
                        step = 1E5)
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
        min = 1,
        max = 3000,
        value = 1000,
        step = 50
      ),
      shiny::sliderInput(
        inputId = 'cofactor.yaxis',
        label = "Cofactor: Y-axis",
        min = 1,
        max = 3000,
        value = 1000,
        step = 50
      )
    )
  })
  cyto.plot1<-shiny::fluidRow(
    shinydashboard::box(
      title=NULL,
      shiny::plotOutput("ggbivariate_plot1"),
      width=8
    )
  )
  # axis.click.select<-shiny::fluidRow(
  #   shinydashboard::box(
  #     title="Axis (X,Y) Selection",
  #     plotly::plotlyOutput("plotly_heat",height="150px"),
  #     width=8
  #   )
  # )
  ##
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = paste("Cyto Plot"),
                                    disable = F),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        marker.menu,
        shinydashboard::menuItemOutput('asinh.menu')
      )
    ),
    shinydashboard::dashboardBody(
      cyto.plot1,
      axis.click.select
    )
  )

  server <- function(input, output) {
    ##
    if(asinh.view==TRUE){
      output$asinh.menu <- asinh.menu
    }
    ##
    row.index <- shiny::reactive({
      set.seed(1337)
      i <- sample(1:total.rows,input$rowsamp)
      return(i)
    })
    ##
    ggbivariate_plot1 <- shiny::reactive({
      if(asinh.view==FALSE){
        p.tmp <- gg.func.bivariate(dat[row.index(),],
                                   x = !!ggplot2::sym(input$marker1),
                                   y = !!ggplot2::sym(input$marker2))
      }else if(asinh.view==TRUE){
        shiny::req(input$asinh.applied)
        if(input$asinh.applied=='Yes'){
          p.tmp <- gg.func.bivariate(dat[row.index(),],
                                     x = asinh(!!ggplot2::sym(input$marker1)/input$cofactor.xaxis),
                                     y = asinh(!!ggplot2::sym(input$marker2)/input$cofactor.yaxis)
          )
        }else if(input$asinh.applied=='No'){
          p.tmp <- gg.func.bivariate(dat[row.index(),],
                                     x = !!ggplot2::sym(input$marker1),
                                     y = !!ggplot2::sym(input$marker2)
          )
        }
      }
      p.tmp +
        ggplot2::labs(title = "All Events",
                      subtitle = paste(length(row.index()), "of", total.rows, "displayed")) +
        ggplot2::xlab(input$marker1) +
        ggplot2::ylab(input$marker2)
    })
    ##
    output$ggbivariate_plot1 <- shiny::renderPlot({
      ggbivariate_plot1()
    })
    ##
    output$plotly_heat <- plotly::renderPlotly(p)
    ##
    clicks <- shiny::reactiveValues(dat = data.frame(marker1 = NA, marker2 = NA))
    click <- shiny::reactive({
      plotly::event_data("plotly_click", priority = 'event', source = 'axis.selection')
    })
    #
    shiny::observeEvent(eventExpr = click(),{
      if(is.na(clicks$dat$marker1)&is.na(clicks$dat$marker2)){
        clicks$dat$marker1 <- click()$x
      }else if(is.na(clicks$dat$marker2)){
        clicks$dat$marker2 <- click()$x
        shiny::updateSelectInput(inputId = 'marker1',
                                 selected = clicks$dat$marker1)
        shiny::updateSelectInput(inputId = 'marker2',
                                 selected = clicks$dat$marker2)
        clicks$dat$marker1 <- NA
        clicks$dat$marker2 <- NA
      }
    })
  }
  ##
  shiny::shinyApp(ui, server)
}

gg.func.bivariate <- function(dat,...,bins=100,fill.limits=c(0,50)){
  ggplot2::ggplot(dat,ggplot2::aes(...)) +
    ggplot2::geom_hex(bins = bins) +
    viridis::scale_fill_viridis(option = "plasma", limits = fill.limits, oob = scales::squish)
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

generate_cluster_medians<-function(dt,use.scale.func=T){
  if(!data.table::is.data.table(dt)&'cluster' %in% names(dt)){
    stop("Need data.table with a 'cluster' column")
  }
  cluster.medians <- dt[ , .SD, .SDcols = is.numeric][, lapply(.SD, stats::median), keyby = cluster][,-'cluster',with=F]
  if(use.scale.func){
    cluster.medians[,names(cluster.medians) := lapply(.SD, function(x) (x - mean(x))/stats::sd(x))]
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
