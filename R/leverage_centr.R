#' Calculate leverage centrality for all nodes for multiple subjects.
#'
#' This function calculates the leverage centrality for a single graph.
#' @param graphs An undirected igraph object.
#' @param weighted By default weighted=FALSE, but can also be set to weighted=TRUE.
#' @return A matrix of the eigenvector centralities of each node for each subject.
#' @export
#' @author Alex Upton, Brandon Vaughan
#' @details  Leverage centrality is a new measure of node centrality put forth by Joyce et al (2010). Developed specifically for use in functional brain network analysis, this measure has a favorable property of attempting to chracterize information flow in an undirected functional brain network. Leverage centrality compares the degree of a node to the degree of all its neighbors. A high degree node is not highly central according to leverage if all of its neighbors are also high degree. Furthermore, leverage centrality does not assume information flows along the shortest path or in a serial fashion, but rather focuses on the disparity in node degrees in a small neighborhood to quantify consolidation and dissemination of information locally (Joyce et al, 2010). A node with a high degree/strength is not necessarily one with a high centrality value. Leverage centrality defines centrality as having high degree/strength relative to the degree/strength of a node's neighbors. If the neighbors are also of high degree/strength the node is not considered a central node (hub). Keeping with its biological inspiration, leverage centrality does not assume information in a network flows in a serial fashion or only along the shortest path, but on how information flows within a local neighborhood of nodes.
#'
#' The formula for leverage centrality is given below, where k_i is the degree/strength of a node and k_j is the degree/strength of its neighbors.
#'
#' \eqn{L = \frac{1}{k_i} \sum_N_i \frac{k_i-k_j}{k_i+k_j}}
#'
#' For more information on the mathematics of leverage centrality see Vargas et al (2017). This function was originally written by Alex Upton and can be found on the igraph wiki (see references). What this function changes is replacing NaN with zero in the case where a node without any connections after thresholding results in divison by zero. This function also expands leverage centrality to weighted undirected networks. If doing this, be sure to threshold the correlation matrix per Joyce et al (2010).
#'
#' @examples
#' leverage_centr(binary.graphs, weighted=FALSE)
#' leverage_centr(weighted.graphs, weighted=TRUE)
#'
#' @references
#' Joyce, K. E., Laurienti, P. J., Burdette, J. H., & Hayasaka, S. (2010). A New Measure of Centrality for Brain Networks. PLoS ONE, 5(8). doi:10.1371/journal.pone.0012200
#'
#' Vargas, R., Waldron, A., Sharma, A., Flórez, R., & Narayan, D. A. (2017). A graph theoretic analysis of leverage centrality. AKCE International Journal of Graphs and Combinatorics, 14(3), 295-306. doi:10.1016/j.akcej.2017.05.001
#'
#' http://igraph.wikidot.com/r-recipes#toc10

leverage_centr = function(graphs, weighted=FALSE){
  lev = pbapply::pbsapply(graphs, function(g) leverage_centr_single(g, weighted=weighted))
  lev = t(lev)
  return(lev)
}
