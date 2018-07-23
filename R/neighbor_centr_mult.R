#' Calculate neighbor centrality for multiple graphs
#'
#' This calculates the average strength or degree of each node's connections for multiple networks.
#' @param graph A list of networks in matrix format or as an igraph object. Can be weighted (and signed) or binary.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @return A matrix of the average strength/degree of a node's neighbors for each graph.
#' @export
#' @author Brandon Vaughan
#'
#' @details This is an alternative to igraph's knn function. The name was changed from knn to avoid
#' confusion with functions of the same name relating to the machine learning method k-nearest neighbors.
#' The neighbor centrality of node is the mean strength or degree of all of its neighbors (connections).
#'
#' @examples
#' neighbor_centr_mult(graphs)
#'
#' @seealso
#' \code{\link[rsfcNet]{strength_multiple}}
#' \code{\link[rsfcNet]{degree_centr}}
#' \code{\link[rsfcNet]{neighbor_centr}}
#' \code{\link[igraph]{knn}}
#'
#'
#' @references
#'
#' Barrat, A., Barthelemy, M., Pastor-Satorras, R., Vespignani, A. (2004). The architecture of complex weighted networks, Proc. Natl. Acad. Sci. USA 101, 3747
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Node Degree and Strength. Chapter 4. Fundamentals of Brain Network Analysis, 115-136. doi:10.1016/B978-0-12-407908-3.00004-2
#'
#' Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. NeuroImage, 56(4), 2068-2079. doi:10.1016/j.neuroimage.2011.03.069
#`
neighbor_centr_mult = function(graphs, col.names=NULL, row.names=NULL) {
  nn_centr = pbapply::pbsapply(graphs, function(G) neighbor_centr(G))
  nn_centr = t(nn_centr)
  colnames(nn_centr) = col.names
  rownames(nn_centr) = row.names
  return(nn_centr)
}

