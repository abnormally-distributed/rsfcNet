#' Calculate closeness centrality for a list of graphs.
#'
#' This function is a convenience wrapper for igraph's closeness function and takes as an input a list of igraph objects.
#' @param graphs a list of igraph objects.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @return A matrix of the closeness centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#'
#' @details Closeness centrality gives the average shortest path length
#' betwen a node i and all other nodes.
#'
#' \eqn{\frac{N-1}{\sum_{j\neq i}^N d(n_i,e_{ij})}}
#'
#' @examples
#' closeness = closeness_centr_mult(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @seealso
#'
#' \code{\link[rsfcNet]{closeness_centr_mult}}
#' \code{\link[rsfcNet]{curent_centr}}
#'
#' @references
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Freeman, L.C. (1979). Centrality in Social Networks I: Conceptual Clarification. Social Networks, 1, 215-239.
#'

closeness_centr_mult = function (graphs, col.names = NULL, row.names = NULL)
{
  closeness = pbapply::pbsapply(graphs, function(x) igraph::closeness(x,normalized = FALSE, weights= abs(E(x)$weight)))
  closeness = t(closeness)
  colnames(closeness) = col.names
  rownames(closeness) = row.names
  closeness = as.matrix(closeness)
  return(closeness)
}


#' Calculate closeness centrality for a single graph.
#'
#' This function is a convenience wrapper for igraph's closeness function for consistency with rsfcNet's naming conventions.
#' @param graph A network as an igraph object or matrix.
#' @param weights a vector of edge weights. Defaults to the absolute value of the input network.
#' @return A vector of centrality scores for each node.
#' @export
#' @author Brandon Vaughan
#'
#' @details Closeness centrality gives the average shortest path length
#' betwen a node i and all other nodes.
#'
#' \eqn{\frac{N-1}{\sum_{j\neq i}^N d(n_i,e_{ij})}}
#'
#' @examples
#' closeness = closeness_centr_mult(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @seealso
#' \code{\link[rsfcNet]{closeness_centr_mult}}
#' \code{\link[rsfcNet]{current_centr}}
#'
#'
#' @references
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Freeman, L.C. (1979). Centrality in Social Networks I: Conceptual Clarification. Social Networks, 1, 215-239.
#'
closeness_centr = function(graph, weights = abs(E(graph)$weight)) {
  if (igraph::is.igraph(graph)==FALSE){
    graph = igraph::graph.adjacency(graph, mode="undirected", weighted=T, diag=FALSE)
  }
  igraph::closeness(graph,normalized = FALSE, weights= weights)
}
