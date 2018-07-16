#' Calculate closeness centrality for all nodes for multiple subjects.
#'
#' This function is a convenience wrapper for igraph's closeness function and takes as an input a list of igraph objects. ECV is calculated on the absolute values of the correlation matrix per Lohmann et al, 2010. ECV is standardized subject-wise to 0-1 range for consistency.
#' @param graphs a list of igraph objects.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @return A matrix of the closeness centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#'
#' @details Closeness centrality (in the normalized form given here) gives the average shortest path length
#' betwen a node i and all other nodes.
#'
#' \eqn{\frac{N-1}{\sum_{j\neq i}^N d(n_i,e_{ij})}}
#'
#' @examples
#' closeness = closeness_centr(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @references
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Freeman, L.C. (1979). Centrality in Social Networks I: Conceptual Clarification. Social Networks, 1, 215-239.
#'

closeness_centr = function (graphs, col.names = NULL, row.names = NULL)
{
  closeness = pbapply::pbsapply(graphs, function(x) igraph::closeness(x,normalized = scale, weights= abs(E(x)$weight)))
  closeness = t(closeness)
  colnames(closeness) = col.names
  rownames(closeness) = row.names
  closeness = as.matrix(closeness)
  return(closeness)
}
