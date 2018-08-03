#' Calculate the global average clustering coefficient for a single graph
#'
#' Calculates the average clustering coefficient for a single graph using any of the variants
#' of weighted, signed, clustering coefficient.
#'
#' @param graph A graph as an igraph object or matrix..
#' @param method Either "Constantini" (the default) or "Zhang"
#'
#'
#' @details One straightofrward measure of the global clustering coefficient was put forth by Watts and Strogatz
#' (1998) as simply the average clustering coefficient across all nodes in a network.
#' @export
#' @return The global clustering coefficient for a graph
#' @author Brandon vaughan
#'
#' @seealso
#' \link[rsfcNet]{local_trans}
#' \link[igraph]{transitivity}
#' \link[rsfcNet]{clustcoef_signed}
#'
#' @references
#' Watts, D. J.  and Strogatz, S. (1998). Collective dynamics of 'small-world' networks. Nature. 393 (6684): 440–442. doi:10.1038/30918
#'
#'

global_clust = function(graph, method="Constantini") {
  c = clustcoef_signed(graph, method=method)
  mean(c)
}


#' Calculate the global average clustering coefficient for a list of graphs.
#'
#' Calculates the average clustering coefficient for a list of graphs using any of the variants
#' of weighted, signed, clustering coefficient.
#'
#' @param graphs A list of igraph objects or connectivity matrices.
#' @param method Either "Constantini" (the default) or "Zhang"
#'
#' @details One straightofrward measure of the global clustering coefficient was put forth by Watts and Strogatz
#' (1998) as simply the average clustering coefficient across all nodes in a network.
#' @export
#' @return A vector of global clustering coefficients for the graphs
#' @author Brandon vaughan
#'
#' @seealso
#' \link[rsfcNet]{global_trans}
#' \link[rsfcNet]{global_clust}
#' \link[rsfcNet]{clustcoef_signed_mult}

global_clust_mult = function(graphs,method= "Constantini") {
  c = clustcoef_signed_mult(graphs, method=method)
  rowMeans(c)
}
