#' Calculate degree centrality for all nodes for multiple subjects.
#'
#' This function is a convenience wrapper for igraph's degree function and takes as an input a list of igraph objects.
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
#' \code{\link[rsfcNet]{strength_multiple}}
#' \code{\link[rsfcNet]{strength_signed}}
#' \code{\link[igraph]{strength}}
#'
#' @examples
#' degree = degree_centr(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @references
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Node Degree and Strength. Chapter 4. Fundamentals of Brain Network Analysis, 115-136. doi:10.1016/B978-0-12-407908-3.00004-2
#'
degree_centr = function (graphs, col.names = NULL, row.names = NULL)
{
  degrees = pbapply::pbsapply(graphs, function(x) igraph::degree(x,normalized = FALSE))
  degrees = t(degree)
  colnames(degrees) = col.names
  rownames(degrees) = row.names
  degree = as.matrix(degrees)
  return(degrees)
}
