#' @title Interactive asinh transformation of cytometry data for determining optimal cofactor values
#'
#' @param dt `data.table` containing raw, un-transformed cytometry data; as returned from `fcs.to.dt`
#' @param cofactors.file.path if provided, will save (.RDS) the cofactor values to the specified \code{file.path}
#' @param filename.suffix Character string; if provided, will be appended to the end of the file name
#'
#' @return \code{shinyapp}; if \code{cofactors.file.path} is provided, a \code{data.table} of cofactor values will be saved for use in transforming raw data values.
#' @export
#'
cofactor.viewer<-function(dt,cofactors.file.path=NULL,filename.suffix=NULL){
  ##
  col.classes<-unlist(sapply(dt,class))
  col.vars<-sapply(unique(col.classes),function(i) names(col.classes[col.classes %in% i]),simplify = F)
  fluors <- grep("Time|FSC|SSC|node|cluster", col.vars$numeric, value = T, invert = T)
  cuts<-Reduce(union,dt[, lapply(.SD, function(x){q<-stats::quantile(x,probs=c(0.0001,.9999));which(x<q[1]|x>q[2])}), .SDcols = fluors])
  if('sample.id' %in% col.vars$character){
    samples<-c('combined',sort(unique(dt[['sample.id']])))
  }else{
    message("Sample selection requires a 'sample.id' column (character).")
  }
  ##
  slider.vals <- stats::setNames(nm = c("min","max","value","step"), c(100, 10000, 1000, 50))
  ##
  p1 <- shinydashboard::box(
    collapsible = F, title = "Asinh Transformed Data",
    shiny::plotOutput("ggbivariate_plot1"),
    width = 6
  )
  ##
  p2 <- shinydashboard::box(
    collapsible = F, title = "Markers / Cofactors",
    shiny::selectInput(inputId = "sample.id",
                       label = "Sample:",
                       choices = samples, selected = samples[1]),
    shiny::selectInput(inputId = "marker1",
                       label = "Marker (x):",
                       choices = col.vars$numeric, selected = col.vars$numeric[1]),
    shiny::sliderInput(inputId = "cofactor.x",
                       label = "Cofactor: X-axis",
                       min = slider.vals[["min"]], max = slider.vals[["max"]],
                       value = slider.vals[["value"]], step = slider.vals[["step"]]),
    shiny::selectInput(inputId = "marker2",
                       label = "Marker (y):",
                       choices = col.vars$numeric, selected = col.vars$numeric[2]),
    shiny::sliderInput(inputId = "cofactor.y",
                       label = "Cofactor: Y-axis",
                       min = slider.vals[["min"]], max = slider.vals[["max"]],
                       value = slider.vals[["value"]], step = slider.vals[["step"]]),
    shiny::actionButton("save.cofactors","Save Cofactors"),
    # shiny::actionButton("stop","Quit App -- Return Cofactors"),
    width = 3
  )
  ##
  p3 <- shinydashboard::box(
    collapsible = F, title = "Cofactors (Saved):",style='height:800px;overflow-y: scroll',
    shiny::tableOutput('cofactors'),
    width = 3)
  ##
  header<-shinydashboard::dashboardHeader(title="Asinh Transformation")
  sidebar<-shinydashboard::dashboardSidebar(collapsed = T)
  body<-shinydashboard::dashboardBody(shiny::fluidRow(p1,p2,p3))
  ##
  ui <- shinydashboard::dashboardPage(
    header,
    sidebar,
    body
  )
  ##
  server <- function(input, output,session) {
    ##
    ggbivariate_plot1 <- shiny::reactive({
      p <- gg.func.bivariate(
        dat=if(input$sample.id=='combined'){
          dt[-cuts]
        }else{
          dt[-cuts][get("sample.id") == input$sample.id]
        },
        x = asinh(!!ggplot2::sym(input$marker1)/input$cofactor.x),
        y = asinh(!!ggplot2::sym(input$marker2)/input$cofactor.y)
      )
      p<-p+ggplot2::xlab(input$marker1) + ggplot2::ylab(input$marker2)
      return(p)
    })
    output$ggbivariate_plot1 <- shiny::renderPlot({
      ggbivariate_plot1()
    })
    ##
    rv <- shiny::reactiveValues(cofactors=data.frame(marker=fluors, cofactor=slider.vals[["value"]]))

    shiny::observeEvent(input$save.cofactors, {
      rv$cofactors[rv$cofactors$marker==input$marker1,'cofactor'] <- input$cofactor.x
      rv$cofactors[rv$cofactors$marker==input$marker2,'cofactor'] <- input$cofactor.y
    })
    ##
    output$cofactors <- shiny::renderTable(rv$cofactors,striped = T,digits = 0)
    ##
    # shiny::observeEvent(input$stop, {
    #   shiny::stopApp()
    # })
    ##
    if(!is.null(cofactors.file.path)){
      session$onSessionEnded(function(){
        shiny::isolate(
          saveRDS(data.table::data.table(rv$cofactors),
                  file = file.path(cofactors.file.path,
                                   ifelse(is.null(filename.suffix),
                                          sprintf("%s_%s.RDS","cofactors",Sys.Date()),
                                          sprintf("%s_%s_%s.RDS","cofactors",filename.suffix,Sys.Date())
                                   )
                  )
          )
        )
      })
    }
    ##
  }
  ##
  shiny::shinyApp(ui, server)
}
