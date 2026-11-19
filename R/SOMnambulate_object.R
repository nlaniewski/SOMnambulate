#' @title Build a Self-organizing Map (SOM)
#'
#' @param flowstate An object of class flowstate; as returned from `[flowstate:::flowstate_object]`.
#' @param dims Character vector; dimensions/variables/columns -- input features.
#' @param seed.val An elite random seed value -- `1337` -- used to control otherwise random starts.
#' @param ... Further arguments passed to [flowsomlite::SOM].
#'
#' @returns A list containing [flowsomlite::SOM]-specific parameters/results. Of primary importance is the list element `[['codes']]`.
#' @export
#'
#'
som <- function(
    flowstate,
    dims,
    seed.val = 1337,
    ...
){
  ##
  flowstate_name <- deparse(substitute(flowstate))
  if(!data.table::is.data.table(flowstate)){
    stop(sprintf("'%s' must be a data.table", flowstate_name), call. = F)
  }
  if(!all(dims %in% names(flowstate))){
    dims.not <- dims[!dims %in% names(flowstate)]
    dims.not <- paste0(dims.not, collapse = " ; ")
    stop(sprintf("%s ('dims')  not in '%s'", dims.not, flowstate_name), call. = F)
  }
  ##
  som <- flowsomlite::SOM(
    data = flowstate[, as.matrix(.SD), .SDcols = dims],
    seed.val = seed.val,
    ...
  )
  ##
  som[['codes.dt']] <- data.table::data.table(
    som$codes,
    node = seq_len(som$nNodes)
  )
  ##
  som[['codes.dt']][
    ,
    j = c('umap.1', 'umap.2') := {
      uwot::umap(
        X = apply(som$codes, 2, function(x) (x - min(x)) / diff(range(x))),
        seed = seed.val
      ) |> as.data.frame() |> as.list()
    }
  ]
  ##
  return(som)
}

#' @title Cluster a Self-organizing Map (SOM)
#'
#' @param som The return of [som].
#' @param k.name Character string; a suffix appended to the otherwise generically named `node` and `cluster` columns.
#' @param k Integer; desired number of clusters.
#' @param plot Logical -- default `FALSE`; plot the [stats::hclust] result.
#'
#' @returns UPDATES BY REFERENCE (invisibly returns `som`); the features in `som[['codes']]` are clustered, the (factored) result of which is added in-place to `som[['codes.dt']]`.
#' @export
#'
#'
som_k <- function(
    som,
    k.name,
    k,
    plot = F
)
{
  .dist <- stats::dist(som$codes, method = "euclidean")
  hc <- stats::hclust(.dist, method = "ward.D2")

  if(plot){
    plot(hc)
    stats::rect.hclust(hc, k = k, border = "red")
  }

  clusters <- factor(stats::cutree(hc, k = k))

  k.name <- sprintf("cluster_%s", k.name)

  som[['codes.dt']][
    ,
    (k.name) := clusters
  ]

  invisible(som)
}

#' @title Map individual data points to a SOM.
#'
#' @param som The return of [som].
#' @param flowstate An object of class flowstate; as returned from `[flowstate:::flowstate_object]`.
#'
#' @returns UPDATES BY REFERENCE (invisible returns `flowstate`); `node`, `node_dist`, and `cluster` columns are added in-place to `flowstate`.
#' @export
#'
#'
som_map <- function(
    som,
    flowstate
){
  flowstate[
    ,
    j = c('node', 'node_dist') := {
      flowsomlite::MapDataToCodes(
        codes = som[['codes']],
        newdata = as.matrix(.SD)
      ) |> as.data.frame() |> as.list()
    },
    .SDcols = colnames(som[['codes']])
  ]

  invisible(flowstate)
}
