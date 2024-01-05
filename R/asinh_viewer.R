#' @title Interactive asinh transformation of cytometry data
#'
#' @param dat \code{data.table} containing raw, untransformed cytometry data
#'
#' @return \code{shinyapp}; \code{asinh.cofactors} - a named vector of cofactor values for use in transforming data - will be assigned to the environment if the \code{actionButton} is clicked.
#' @export
#'
cofactor.viewer<-function(dat){
  ##
  vars<-stats::setNames(nm=c("num","fac","char"),c("is.numeric","is.factor","is.character"))
  cols<-sapply(vars,function(i) stringr::str_sort(names(dat[,.SD,.SDcols = get(i)]),numeric = T))
  fluors<-grep("Time|FSC|SSC",cols$num,value = T,invert = T)
  cuts<-Reduce(union,dat[, lapply(.SD, function(x){q<-stats::quantile(x,probs=c(0.0001,.9999));which(x<q[1]|x>q[2])}), .SDcols = fluors])
  if('sample.id' %in% cols$char){
    samples<-c('combined',sort(unique(dat[['sample.id']])))
  }else{
    message("Sample selction requires a 'sample.id' column (character strings).")
  }
  ##
  slider.vals <- stats::setNames(nm = c("min","max","value","step"), c(100, 10000, 1000, 50))
  ##
  p1 <- shinydashboard::box(
    collapsible = F, title = "Asinh Transformed Data",
    shiny::plotOutput("ggbivariate_plot1"),
    width = 8
  )
  ##
  p2 <- shinydashboard::box(
    collapsible = F, title = "Markers / Cofactors",
    shiny::selectInput(inputId = "sample.id",
                       label = "Sample:",
                       choices = samples, selected = samples[1]),
    shiny::selectInput(inputId = "marker1",
                       label = "Marker (x):",
                       choices = cols$num, selected = cols$num[1]),
    shiny::sliderInput(inputId = "cofactor.x",
                       label = "Cofactor: X-axis",
                       min = slider.vals[["min"]], max = slider.vals[["max"]],
                       value = slider.vals[["value"]], step = slider.vals[["step"]]),
    shiny::selectInput(inputId = "marker2",
                       label = "Marker (y):",
                       choices = cols$num, selected = cols$num[2]),
    shiny::sliderInput(inputId = "cofactor.y",
                       label = "Cofactor: Y-axis",
                       min = slider.vals[["min"]], max = slider.vals[["max"]],
                       value = slider.vals[["value"]], step = slider.vals[["step"]]),
    shiny::actionButton("save.cofactors","Save Cofactors"),
    width = 4
  )
  ##...
  header<-shinydashboard::dashboardHeader(title="Asinh Transformation")
  sidebar<-shinydashboard::dashboardSidebar(collapsed = T)
  body<-shinydashboard::dashboardBody(shiny::fluidRow(p1,p2))
  ##
  ui <- shinydashboard::dashboardPage(
    header,
    sidebar,
    body
  )
  ##
  server <- function(input, output) {
    ##
    ggbivariate_plot1 <- shiny::reactive({
      p <- gg.func.bivariate(
        dat=if(input$sample.id=='combined'){
          dat[-cuts]
        }else{
          dat[-cuts][get("sample.id") == input$sample.id]
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
    shiny::observeEvent(input$save.cofactors,{
      cfs<-stats::setNames(c(input$cofactor.x,input$cofactor.y),nm=c(input$marker1,input$marker2))
      if(exists('asinh.cofactors')){
        asinh.cofactors<-c(asinh.cofactors,cfs)
      }else{
        asinh.cofactors<-cfs
      }
      asinh.cofactors<<-asinh.cofactors[unique(names(asinh.cofactors))]
    })
    ##
  }
  ##
  shiny::shinyApp(ui, server)
}
