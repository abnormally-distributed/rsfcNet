#' Calculate the average connectivity of each node's neighbors
#'
#' This calculates the average strength or degree of each node's connections for a single network.
#' @param graph A network in matrix format or as an igraph object. Can be weighted (and signed) or binary.
#' @return A vector of the average strength/degree of a node's neighbors.
#' @export
#' @author Brandon Vaughan
#'
#' @details This is an alternative to igraph's knn function. The name was changed from knn to avoid
#' confusion with functions of the same name relating to the machine learning method k-nearest neighbors.
#' The neighbor centrality of node is the mean strength or degree of all of its neighbors (connections).
#'
#'
#' @examples
#' neighbor_centr(graph)
#'
#' @seealso
#' \code{\link[rsfcNet]{strength_signed}}
#' \code{\link[igraph]{degree}}
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
neighbor_centr = function(graph) {
if (is.igraph(graph)=="TRUE"){
  graph = as.matrix(igraph::as_adjacency_matrix(graph, edges = FALSE, attr = "weight", sparse=TRUE)) }
  graph.pos = graph
  graph.neg = graph
  graph.pos[which(graph<0)] <-0
  graph.neg[which(graph>0)] <-0
  graph.neg = graph.adjacency(graph.neg, mode="undirected", weighted=T, diag=TRUE)
  graph.pos= graph.adjacency(graph.pos, mode="undirected", weighted=T, diag=TRUE)
  knn.neg = knn(graph.neg, vids = V(graph.neg), weights = E(graph.neg)$weight)$knn
  knn.neg[which(is.na(knn.neg))] <- 0
  knn.pos[which(is.na(knn.pos))] <- 0
  knn.pos = knn(graph.pos, vids = V(graph.pos), weights = E(graph.pos)$weight)$knn
  knn.neg = (knn.neg/(knn.pos+knn.neg))*knn.neg
  knn.neg[which(is.na(knn.neg))] <- 0
  knn.star = knn.pos-knn.neg
  knn.star[which(is.na(knn.star))] <- 0
  return(knn.star)
}

