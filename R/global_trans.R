#' Calculate the global transitivity of a list of graphs
#'
#' This function is a convenience wrapper for the transitivity function in igraph. It is used to calculate the global clustering coefficient of a list of unweighted graphs or a thresholded weighted graph
#' @param graphs a list of igraph objects.
#' @return A vector of global transitivity scores for each subject
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' Global transitivity is the ratio of complete triangles to all triplets (three nodes connected by at least two triplets, ie a "V" shape or a complete triangle) in the network. It is a variant of the global clustering coefficient (Wasserman & Faust, 1994).
#'
#' @examples
#' GTR = global_trans(graphs)
#'
#' @references
#' Wasserman, S., and Faust, K. (1994). Social Network Analysis: Methods and Applications. Cambridge: Cambridge University Press.

global_trans = function(graphs) {
pbapply::pbsapply(graphs, function(g) igraph::transitivity(g, type="global", isolates="zero"))
}

