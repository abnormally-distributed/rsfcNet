#' Calculate the local transitivity of a list of graphs
#'
#' This function is a convenience wrapper for the transitivity function in igraph. It is used to calculate the local clustering coefficient for each node in a list of unweighted graphs or a thresholded weighted graph.
#' @param graphs a list of igraph objects. Can be weighted or unweighted, but will not perform well without thresholding applied to the weighted graph.
#' @return A matrix oflocal transitivity scores for all nodes in each subject
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' The local transitivity of a node is the fraction of edges a node forms with its neighbors out of the the number of edges it would take to make complete triangles. Barratt's formula is used, which gives the normal local transitivity if the graph is unweighted but also allows for weighted edges (Barratt et al, 2004).
#'
#' @examples
#' LTR = local_trans(graphs)
#'
#' @references
#' Barrat A, Barthélemy M, Pastor-Satorras R, Vespignani A. The architecture of complex weighted networks. Proceedings of the National Academy of Sciences of the United States of America. 2004;101(11):3747-3752. doi:10.1073/pnas.0400087101.
#'
> local_trans = function (graphs) 
{
     pbapply::pbsapply(graphs, function(g) igraph::transitivity(g, type = "barrat", isolates = "zero"))
 }

