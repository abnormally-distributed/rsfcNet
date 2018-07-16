#' Calculate leverage centrality for all nodes for multiple subjects.
#'
#' This function calculates the leverage centrality for a single graph.
#' @param graphs a list of igraph objects.
#' @param weighted By default weighted=FALSE, but can also be set to weighted=TRUE.
#' @param parallel Should multiple cores be used? Defaults to FALSE. If TRUE, progress bar is not displayed. This is normal.
#' @param cores How many cores should be used? Defaults to recommended 1 less than number of CPU cores.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @return A matrix of the leverage centralities of each node for each subject.
#' @export
#' @author Alex Upton, Brandon Vaughan
#' @details  Leverage centrality is a new measure of node centrality put forth by Joyce et al (2010). Developed specifically for use in functional brain network analysis, this measure has a favorable property of attempting to chracterize information flow in an undirected functional brain network. Leverage centrality compares the degree of a node to the degree of all its neighbors. A node with a high degree/strength is not necessarily one with a high centrality value. Leverage centrality defines centrality as having high degree/strength relative to the degree/strength of a node's neighbors. If the neighbors are also of high degree/strength the node is not considered a central node. Keeping with its biological inspiration, leverage centrality does not assume information in a network flows in a serial fashion or only along the shortest path, but on how information flows within a local neighborhood of nodes (Joyce et al, 2010).
#'
#' The formula for leverage centrality is given below, where k_i is the degree/strength of a node and k_j is the degree/strength of its neighbors.
#'
#' \eqn{L = \frac{1}{k_i} \sum_N_i \frac{k_i-k_j}{k_i+k_j}}
#'
#' For more information on the mathematics of leverage centrality see Vargas et al (2017). This function was originally written by Alex Upton and can be found on the igraph wiki (see references). What this function changes is replacing NaN with zero in the case where a node without any connections after thresholding results in divison by zero. This function also expands leverage centrality to weighted undirected networks. If doing this, be sure to threshold the correlation matrix per Joyce et al (2010).
#'
#' @examples
#' leverage_centr_mult(binary.graphs, weighted=FALSE)
#' leverage_centr_mult(weighted.graphs, weighted=TRUE)
#'
#' @seealso
#' \code{\link[rsfcNet]{leverage_centr}}
#' \code{\link[rsfcNet]{laplace_centr}}
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#'
#' @references
#' Joyce, K. E., Laurienti, P. J., Burdette, J. H., & Hayasaka, S. (2010). A New Measure of Centrality for Brain Networks. PLoS ONE, 5(8). doi:10.1371/journal.pone.0012200
#'
#' Vargas, R., Waldron, A., Sharma, A., Flórez, R., & Narayan, D. A. (2017). A graph theoretic analysis of leverage centrality. AKCE International Journal of Graphs and Combinatorics, 14(3), 295-306. doi:10.1016/j.akcej.2017.05.001
#'
#' http://igraph.wikidot.com/r-recipes#toc10
#'
leverage_centr_mult = function(graphs, weighted=FALSE, col.names=NULL, row.names=NULL, parallel=FALSE, cores=NA){

  if (parallel==FALSE) {
    lev = pbapply::pbsapply(graphs, function(g) leverage_centr(g, weighted=weighted))
    rownames(lev) = row.names
    colnames(lev) = col.names
    return(lev)
  } else if (parallel==TRUE) {

    if (is.na(cores)==TRUE) {
      cl = as.integer(parallel::detectCores()-1)
      cl = makeCluster(cl)
    } else if (is.na(cores)==FALSE) {
      cl = as.integer(cores)
      cl = makeCluster(cl)
    }

    clusterEvalQ(cl, {
      library(igraph)
    })

    lev = parallel::parSapply(cl,graphs, function(x) rsfcNet::leverage_centr(x, weighted=weighted))
    stopCluster(cl)
    colnames(lev) = col.names
    rownames(lev) = row.names
    lev = as.matrix(lev)
    return(lev)
  }
}

