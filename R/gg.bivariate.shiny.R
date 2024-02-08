#' @title Create a clickable `plotly::plot_ly` heatmap
#'
#' @param vars Character string ('variables') of class `numeric` or `factor`
#' @param type `xy` returns a clickable heatmap of numeric variables; `factor` returns a clickable heatmap of factored variables; `cluster` (if column is present) returns a clickable heatmap of cluster numbers.
#'
#' @return a `plotly::plot_ly` heatmap; allows for clickable selection of reactive variables when used with `shiny`
#'
#'
#'
selection.heatmap<-function(vars,type=c('xy','factor','cluster')){
  if(type=='xy'){
    bar1<-rep(c(0,1),length(vars))[seq(vars)]
    if(length(bar1) %%2 == 1){
      bar2<-rep(c(1,0),length(vars))[seq(vars)]
    }else{
      bar2<-rev(bar1)
    }
    bars<-matrix(c(bar1,bar2),ncol=2,dimnames=list(vars,c('click.x','click.y')))
  }else{
    vars<-c('NULL',vars)
    bar1<-rep(c(0,1),length(vars))[seq(vars)]
    bars<-matrix(bar1,ncol=1,dimnames=list(vars,paste0('click.',type)))
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
    title=ifelse(type=='xy',"Axis Selection",ifelse(type=='factor',"Factor Selection","Cluster Selection")),
    xaxis = list(tickfont = list(size = 10),type = "category"),
    yaxis = list(tickfont = list(size = 10),type = "category")
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
          shiny::column(
            width = 4,
            plotly::plotlyOutput(outputId = "axis.select",height = "80vh")
          ),
          shiny::column(
            width = 4,
            plotly::plotlyOutput(outputId = "factor.select",height = "20vh")
          ),
          if(!is.null(cluster.var)){
            shiny::column(
              width= 4,
              plotly::plotlyOutput(outputId = "cluster.select",height = "80vh")
            )
          }
        ),
        width=6
      ),
      shiny::mainPanel(
        shiny::plotOutput(outputId = 'gg',height = "800px"),
        width=6
      )
    )
  )
  #
  server<-function(input, output, session) {
    output$axis.select <- plotly::renderPlotly(selection.heatmap(vars=xy.vars,type='xy'))
    clicks.axis <- shiny::reactiveValues(x=xy.vars[1],y=xy.vars[2])
    clicked.axis <- shiny::reactive({
      plotly::event_data(
        event = "plotly_click",
        source = "xy.select",
        priority = "event")
    })
    shiny::observeEvent(eventExpr = clicked.axis(),{
      if(clicked.axis()$x=='click.x'){
        clicks.axis$x<-clicked.axis()$y
      }else if(clicked.axis()$x=='click.y'){
        clicks.axis$y<-clicked.axis()$y
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
    if(!is.null(cluster.var)){
      output$cluster.select <- plotly::renderPlotly(selection.heatmap(vars=cluster.vals,type='cluster'))
      clicks.cluster <- shiny::reactiveVal('NULL')
      clicked.cluster <- shiny::reactive({
        plotly::event_data(
          event = "plotly_click",
          source = "cluster.select",
          priority = "event")
      })
      shiny::observeEvent(eventExpr = clicked.cluster(),{
        clicks.cluster(clicked.cluster()$y)
      })
    }
    #
    output$gg <- shiny::renderPlot({
      ggplot2::ggplot(
        if(is.null(cluster.var)){
          dt
        }else if(!is.null(cluster.var)&clicks.cluster()=='NULL'){
          dt
        }else if(!is.null(cluster.var)&clicks.cluster()!='NULL'){
          dt[get(cluster.var) %in% clicks.cluster()]
        },
        ggplot2::aes(!!as.name(clicks.axis$x),!!as.name(clicks.axis$y))) +
        #programmatically determine optimal bin number?
        ggplot2::geom_hex(bins=200) +
        viridis::scale_fill_viridis(option='plasma') +
        ggplot2::facet_wrap(ifelse(clicks.factor$x=='none',NA,clicks.factor$x))
    })
  }
  #
  shiny::shinyApp(ui, server)
}
