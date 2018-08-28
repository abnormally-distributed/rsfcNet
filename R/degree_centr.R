#' Calculate degree centrality for a single graph.
#'
#' This function calculates the number of non-zero connections each node in a network has for a single graph.
#' @param graph A network as an igraph object or connectivity matrix.
#' @return A vector of centrality scores for each node.
#' @export
#' @author Brandon Vaughan
#'
#' @details Degree centrality is simply the sum of present connections of node i to each other jth node:
#'
#' \eqn{D_i= \sum_{j \neq i}^N e_{ij}}
#'
#' For the weighted version of degree, see the functions for calculating strength.
#'
#' @seealso
#' \code{\link[rsfcNet]{degree_mult}}
#' \code{\link[rsfcNet]{strength_signed}}
#' \code{\link[igraph]{strength}}
#'
#' @examples
#'
#' # degree = degree_centr(graph)
#'
#' @references
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Node Degree and Strength. Chapter 4. Fundamentals of Brain Network Analysis, 115-136. doi:10.1016/B978-0-12-407908-3.00004-2
#'
degree_centr = function(graph)
{
  if (is.igraph(graph)=="TRUE"){
    graph = as.matrix(as_adj(graph, type="both", edges=FALSE, sparse=TRUE))
  }
  graph = rsfcNet::binarize(graph)
  degrees = colSums(degrees)
  return(degrees)
}


#' Calculate degree centrality for multiple graphs.
#'
#' This function calculates the number of non-zero connections each node in a network has for a list of graphs.
#' @param graphs a list of igraph objects.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @return A matrix of the degree centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#'
#' @details Degree centrality is simply the sum of present connections of node i to each other jth node:
#'
#' \eqn{D_i= \sum_{j \neq i}^N e_{ij}}
#'
#' For the weighted version of degree, see the functions for calculating strength.
#'
#' @seealso
#' \code{\link[rsfcNet]{degree_centr}}
#' \code{\link[rsfcNet]{strength_multiple}}
#'
#' @examples
#'
#' # degree = degree_mult(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @references
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Node Degree and Strength. Chapter 4. Fundamentals of Brain Network Analysis, 115-136. doi:10.1016/B978-0-12-407908-3.00004-2
#'
degree_mult = function (graphs, col.names = NULL, row.names = NULL)
{
  degrees = pbapply::pbsapply(graphs, function(x) rsfcNet::degree_centr(x))
  degrees = t(degree)
  colnames(degrees) = col.names
  rownames(degrees) = row.names
  degree = as.matrix(degrees)
  return(degrees)
}
