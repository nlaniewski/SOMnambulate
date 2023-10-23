#' @title Generate bivariate plot pairs
#'
#' @param plot.dims Character vector; a string of dimensions/markers from which bivarate pairs will be generated
#'
#' @return a list of plot pairs; useful for \code{ggplot2} loops
#'
#'
#'
plot.pairs<-function(plot.dims){
  pair.1 <- plot.dims[seq(1, length(plot.dims), 2)]
  pair.2 <- plot.dims[seq(2, length(plot.dims), 2)]
  if (length(pair.1) != length(pair.2)) {
    pair.2[length(pair.1)] <- pair.2[length(pair.1) - 1]
  }
  pair.list <- mapply(c, pair.1, pair.2, SIMPLIFY = F, USE.NAMES = F)
}
