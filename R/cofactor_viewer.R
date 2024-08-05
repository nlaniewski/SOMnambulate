#' @title Interactive asinh transformation of cytometry data for determining optimal cofactor values
#'
#' @param dt \link[data.table]{data.table} containing raw, un-transformed cytometry data; as returned from \link{fcs.to.dt}.
#' @param cytometry.type One of either 'spectral', 'conventional', or 'mass'; depending on type, a default cofactor will be set.
#' @param cofactors.file.path if provided, will save (.RDS) the cofactor values to the specified \code{file.path}.
#' @param filename.suffix Character string; if provided, will be appended to the end of the file name.
#'
#' @return \code{shinyapp}; if \code{cofactors.file.path} is provided, a \code{data.table} of cofactor values will be saved for use in transforming raw data values.
#' @export
#' @examples
#' \dontrun{
#' dt<-SOMnambulate:::prepared.examples(example.type='dt')
#'
#' #set and save cofactors
#' cofactor.viewer(dt[,.(CD3_170Er,CD4_145Nd,CD8a_162Dy,TCRgd_155Gd,sample.id)],
#' cytometry.type = 'mass',cofactors.file.path=tempdir(),filename.suffix="ECHO")
#'
#' #load saved cofactors
#' readRDS(list.files(tempdir(),full.names=TRUE,pattern='cofactors'))
#' }
#'
cofactor.viewer<-function(dt,cytometry.type=c('spectral','conventional','mass'),cofactors.file.path=NULL,filename.suffix=NULL){
  ##
  col.classes<-unlist(sapply(dt,class))
  col.vars<-sapply(unique(col.classes),function(i) names(col.classes[col.classes %in% i]),simplify = F)
  vars <- grep("Time|FSC|SSC|cluster|node|barcode", col.vars$numeric, value = T, invert = T)
  cuts<-Reduce(union,dt[, lapply(.SD, function(x){q<-stats::quantile(x,probs=c(0.0001,.9999));which(x<q[1]|x>q[2])}), .SDcols = vars])
  samples<-'combined'
  if(!is.null(col.vars$character)&'sample.id' %in% col.vars$character){
    samples<-c(samples,sort(unique(dt[['sample.id']])))
  }else{
    message("Sample selection requires a 'sample.id' column (character).")
  }
  ##
  if(cytometry.type=='spectral'){
    slider.vals <- stats::setNames(nm = c("min","max","value","step"), c(100, 10000, 5000, 50))
  }else if(cytometry.type=='conventional'){
    slider.vals <- stats::setNames(nm = c("min","max","value","step"), c(100, 10000, 1000, 50))
  }else if(cytometry.type=='mass'){
    slider.vals <- stats::setNames(nm = c("min","max","value","step"), c(1, 10, 5, 1))
  }
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
    shiny::actionButton("set.cofactors","Set Cofactors"),
    if(!is.null(cofactors.file.path)){
      shiny::actionButton("save.cofactors","Save Cofactors")
    },
    width = 3
  )
  ##
  p3 <- shinydashboard::box(
    collapsible = F, title = "Cofactors (Set):",style='height:800px;overflow-y: scroll',
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

    rv <- shiny::reactiveValues(cofactors=stats::setNames(nm=vars,rep(slider.vals[["value"]],length(vars))))

    shiny::observeEvent(input$set.cofactors, {
      rv$cofactors[[input$marker1]] <- input$cofactor.x
      rv$cofactors[[input$marker2]] <- input$cofactor.y
    })

    shiny::observeEvent(input$save.cofactors,{
      shiny::isolate(
        saveRDS(rv$cofactors,
                file=file.path(cofactors.file.path,
                               ifelse(is.null(filename.suffix),
                                      sprintf("%s_%s.RDS","cofactors",Sys.Date()),
                                      sprintf("%s_%s_%s.RDS","cofactors",filename.suffix,Sys.Date())
                               )
                )
        )
      )
    })

    output$cofactors <- shiny::renderTable(
      data.frame(marker=names(rv$cofactors),cofactor=rv$cofactors),
      striped = T,digits = 0
    )
  }
  shiny::shinyApp(ui, server)
}
