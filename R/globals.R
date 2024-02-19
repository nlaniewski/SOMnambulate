.datatable.aware = TRUE
utils::globalVariables('.rs.invokeShinyWindowViewer')
#for R CMD check; data.table vars
utils::globalVariables(
  c(
    '.'
    ,'alias'
    ,'barcode'
    ,'barcode.id'
    ,'barcode_node'
    ,'cell.type'
    ,'channels'
    ,'cluster'
    ,'desc'
    ,'maxRange'
    ,'minRange'
    ,'N'
    ,'name'
    ,'node'
    ,'S'
    ,'total'
    ,'type'
    ,'valley'
    ,'value'
    ,'y'
    ,'yield'
    )
)
