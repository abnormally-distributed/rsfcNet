#' Convert list of connectivity matrices to igraph objects
#'
#' This function converts a list of adjacency or correlation matrices to igraph objects.
#' @param cormats a list of matrices
#' @return A list of igraph objects
#' @export
#' @author Brandon Vaughan
#' @examples
#' graphs = mats_to_graphs(cormats)
#'
mats_as_graphs = function(cormats) {
  pbapply::pblapply(cormats, function(x) graph.adjacency(x, weighted=T, mode="undirected", diag=TRUE))
}
