#' Calculate leverage centrality for a single graph..
#'
#' This function calculates the leverage centrality for a single graph.
#' @param graph An undirected igraph object.
#' @param weighted By default TRUE, but can also be set to FALSE.
#' @param strength.star Should the strength* measure of Rubinov & Sporns (2011) be used? Defaults to FALSE. If set to TRUE, only positive weights are used. Only applicable when weighted = TRUE.
#' @return A matrix of the leverage centralities of each node in a subject.
#' @export
#' @author Alex Upton, Brandon Vaughan
#' @details  Leverage centrality is a new measure of node centrality put forth by Joyce et al (2010). Developed specifically for use in functional brain network analysis, this measure has a favorable property of attempting to chracterize information flow in an undirected functional brain network. Leverage centrality compares the degree of a node to the degree of all its neighbors. A node with a high degree/strength is not necessarily one with a high centrality value. Leverage centrality defines centrality as having high degree/strength relative to the degree/strength of a node's neighbors. If the neighbors are also of high degree/strength the node is not considered a central node. Keeping with its biological inspiration, leverage centrality does not assume information in a network flows in a serial fashion or only along the shortest path, but on how information flows within a local neighborhood of nodes (Joyce et al, 2010).
#'
#' The formula for leverage centrality is given below, where k_i is the degree/strength of a node and k_j is the degree/strength of its neighbors.
#'
#' \eqn{L = \frac{1}{k_i} \sum_N_i \frac{k_i-k_j}{k_i+k_j}}
#'
#' For more information on the mathematics of leverage centrality see Vargas et al (2017). This function was originally written by Alex Upton and can be found on the igraph wiki (see references).
#'
#' @examples
#' leverage_centr(binary.graph, weighted=FALSE)
#' leverage_centr(weighted.graph, weighted=TRUE)
#'
#' @references
#' Joyce, K. E., Laurienti, P. J., Burdette, J. H., & Hayasaka, S. (2010). A New Measure of Centrality for Brain Networks. PLoS ONE, 5(8). doi:10.1371/journal.pone.0012200 \cr
#' \cr
#' Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. NeuroImage, 56(4), 2068-2079. doi:10.1016/j.neuroimage.2011.03.069 \cr
#' \cr
#' Vargas, R., Waldron, A., Sharma, A., Flórez, R., & Narayan, D. A. (2017). A graph theoretic analysis of leverage centrality. AKCE International Journal of Graphs and Combinatorics, 14(3), 295-306. doi:10.1016/j.akcej.2017.05.001 \cr
#' \cr
#' http://igraph.wikidot.com/r-recipes#toc10

leverage_centr = function(graph, weighted=TRUE, strength.star = F){
  if (weighted==TRUE){
    if (strength.star){k <- strength_signed(graph)$strength_star}
    else{k <- strength_signed(graph)$positive_strength}
    k <- (k-min(k))/(max(k)-min(k))
    n <- vcount(graph)
    l = sapply(1:n, function(v) { mean((k[v]-k[neighbors(graph,v)]) / (k[v]+k[neighbors(graph,v)])) })
    #l[which(l=="NaN")] <- 0
    return(l)
  }
  else if (weighted==FALSE){
    k <- degree(graph)
    n <- vcount(graph)
    l = sapply(1:n, function(v) { mean((k[v]-k[neighbors(graph,v)]) / (k[v]+k[neighbors(graph,v)])) })
    #l[which(l=="NaN")] <- 0
    return(l)
  }
}
