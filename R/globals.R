.datatable.aware = TRUE
utils::globalVariables('.rs.invokeShinyWindowViewer')
#for R CMD check; data.table vars
utils::globalVariables(
  c(
    '.'
    ,'alias'
    ,'barcode'
    ,'barcode.id'
    ,'barcode.N'
    ,'barcode_node'
    ,'cell.type'
    ,'censor'
    ,'channels'
    ,'cluster'
    ,'d'
    ,'density.x'
    ,'density.y'
    ,'desc'
    ,'f.path'
    ,'maxRange'
    ,'minRange'
    ,'N'
    ,'N.freq'
    ,'name'
    ,'node'
    ,'node_barcode'
    ,'S'
    ,'sample.id'
    ,'Time'
    ,'total'
    ,'type'
    ,'umap.1'
    ,'umap.2'
    ,'variable'
    ,'valley'
    ,'value'
    ,'value.scaled'
    ,'y'
    ,'yield'
    )
)
