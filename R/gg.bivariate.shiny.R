#' @title Create a clickable `plotly::plot_ly` heatmap
#'
#' @param vars Character string ('variables') of class `numeric` or `factor`
#' @param type `xy` returns a clickable heatmap of numeric variables; `factor` returns a clickable heatmap of factored variables.
#'
#' @return a `plotly::plot_ly` heatmap; allows for clickable selection of reactive variables when used with `shiny`
#'
#'
#'
selection.heatmap<-function(vars,type=c('xy','factor')){
  if(type=='xy'){
    bar1<-rep(c(0,1),length(vars))[seq(vars)]
    if(length(bar1) %%2 == 1){
      bar2<-rep(c(1,0),length(vars))[seq(vars)]
    }else{
      bar2<-rev(bar1)
    }
    bars<-matrix(c(bar1,bar2),ncol=2,dimnames=list(vars,c('click.x','click.y')))
  }else if(type=='factor'){
    vars<-c('none',vars)
    bar1<-rep(c(0,1),length(vars))[seq(vars)]
    bars<-matrix(bar1,ncol=1,dimnames=list(vars,'click.factor'))
  }
  plotly.heatmap <- plotly::plot_ly(
    x = colnames(bars),
    y = rev(rownames(bars)),
    z = bars,
    zauto = F,
    type = "heatmap",
    xgap = 2,
    ygap = 2,
    source = paste0(type,".select")
  )
  plotly.heatmap <- plotly::layout(
    plotly.heatmap,
    title=ifelse(type=='xy',"Axis (x,y) Selection","Factor Selection"),
    xaxis = list(tickfont = list(size = 10),type = "category")
  )
  plotly.heatmap <- plotly::config(
    plotly.heatmap,
    displayModeBar=F
  )
  plotly.heatmap <- plotly::layout(
    plotly.heatmap,
    xaxis = list(showspikes = T,spikedash = "longdash", spikemode = "across", spikecolor = "purple",spikethickness = 1),
    yaxis = list(showspike = T, spikedash = "longdash", spikemode = "across", spikecolor = "orange", spikethickness = 1)
  )
  plotly.heatmap <- plotly::hide_colorbar(plotly.heatmap)
  return(plotly.heatmap)
}
#' @title Interactive `shiny::shinyApp`; clickable selection of `ggplot2::ggplot` parameters
#'
#' @param dt a `data.table` as returned from `fcs.to.dt`
#'
#' @return Runs `shiny::shinyApp`
#' @export
#'
#'
gg.bivariate.shiny<-function(dt){
  xy.vars<-colnames(dt)[sapply(dt,class) %in% "numeric"]
  f.vars<-colnames(dt)[sapply(dt,class) %in% "factor"]
  if(any(grepl('cluster',f.vars))){
    cluster.var<-grep('cluster',f.vars,value = T)
    cluster.vals<-dt[,levels(get(cluster.var))]
  }else{
    cluster.var<-NULL
  }
  #
  ui<-shiny::fluidPage(
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fluidRow(
          plotly::plotlyOutput(outputId = "xy.select",height = "80vh")
        ),
        shiny::fluidRow(
          plotly::plotlyOutput(outputId = "factor.select",height = "20vh")
        ),
        if(!is.null(cluster.var)){
          shiny::fluidRow(
            shiny::selectInput(
              inputId = 'cluster',
              label = 'Cluster #',
              choices = cluster.vals,
              selected = NULL
            )
          )
        },
        width=3
      ),
      shiny::mainPanel(
        shiny::plotOutput(outputId = 'gg',height = "800px"),
        width=9
      )
    )
  )
  #
  server<-function(input, output, session) {
    output$xy.select <- plotly::renderPlotly(selection.heatmap(vars=xy.vars,type='xy'))
    clicks.xy <- shiny::reactiveValues(x=xy.vars[1],y=xy.vars[2])
    clicked.xy <- shiny::reactive({
      plotly::event_data(
        event = "plotly_click",
        source = "xy.select",
        priority = "event")
    })
    shiny::observeEvent(eventExpr = clicked.xy(),{
      if(clicked.xy()$x=='click.x'){
        clicks.xy$x<-clicked.xy()$y
      }else if(clicked.xy()$x=='click.y'){
        clicks.xy$y<-clicked.xy()$y
      }
    })
    #
    output$factor.select <- plotly::renderPlotly(selection.heatmap(vars=f.vars,type='factor'))
    clicks.factor <- shiny::reactiveValues(x='none')
    clicked.factor <- shiny::reactive({
      plotly::event_data(
        event = "plotly_click",
        source = "factor.select",
        priority = "event")
    })
    shiny::observeEvent(eventExpr = clicked.factor(),{
      clicks.factor$x<-clicked.factor()$y
    })
    #
    output$gg <- shiny::renderPlot({
      ggplot2::ggplot(dt, ggplot2::aes(!!as.name(clicks.xy$x),!!as.name(clicks.xy$y))) +
        #programmatically determine optimal bin number?
        ggplot2::geom_hex(bins=200) +
        viridis::scale_fill_viridis(option='plasma') +
        ggplot2::facet_wrap(ifelse(clicks.factor$x=='none',NA,clicks.factor$x))
    })
  }
  #
  shiny::shinyApp(ui, server)
}
