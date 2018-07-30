#' Find the closeness vitality for a single graph.
#'
#' Closeness vitality of a node is the change in the sum of distances between all node pairs when that node is removed from the network.
#' @param graph A network as an igraph object or connectivity matrix.
#' @param method One of "regular" (the default) or "current"
#'
#' @return A numeric vector contaning the vitality scores for each node.
#' @author Brandon Vaughan
#' @export
#'
#' @details
#' Closeness vitality of a node is the change in the sum of distances between all node pairs when that node is removed from the network.
#' Either traditional closeness centrality can be used to define the sum of distances (Wiener Index) or the current closeness
#' centrality metric can be substituted with method="current". This function uses the absolute value of the weights, so if you require
#' only positive edges to be considered you must threshold them out of the network.
#'
#' The formula is fairly straightforward:
#'
#' First compute the Wiener Index for the entire graph, which is the sum of the reciprocals of the closeness
#' centrality for each node:
#'
#'
#' \eqn{W=\sum_{i=n}^N \frac{1}{ClC_n}}
#'
#'
#' Next, for each node, calculate the Wiener index again for a subgraph including all nodes but node_i. The vitality
#' of the node is simply the difference between the total Wiener index and the Wiener index for the subgraph.
#'
#'
#' \eqn{Vitality = W_{total} - W_{subgraph}}
#'
#'
#' The vitality of a node measures the increase in transport cost for the whole network
#' when a node is dropped. Negative values indicate a node can be excluded and decrease costs.
#' As the name suggests, the metric measures who in a network is "vital to the operation" and
#' who is an expense to the network (Brandes & Erlebach, 2005). Also see \code{\link[rsfcNet:laplace_centr_mult]{laplacian centrality}}
#' for another measure of node centrality based on the impact to the network after deletion.
#'
#' @references
#' Brandes, U. & Erlebach, T. (2005). Network Analysis: Methodological Foundations, U.S. Government Printing Office.
#'
#' Brandes U., Fleischer D. (2005) Centrality Measures Based on Current Flow. In: Diekert V., Durand B. (eds) STACS 2005. STACS 2005. Lecture Notes in Computer Science, vol 3404. Springer, Berlin, Heidelberg
#'
#' @seealso
#' \code{\link[rsfcNet]{closeness_centr}}
#' \code{\link[rsfcNet]{current_centr}}
#' \code{\link[rsfcNet]{laplace_centr}}
#'
#' @examples
#'
#' **## Not run:**
#' vitality(graph)
#' ## End(**Not run**)
#'
vitality <- function(graph, method="regular") {

  if (method=="current") {
    w <- sum(1/rsfcNet::current_centr(graph))
    count = igraph::vcount(graph)
    nodes <- sapply(1:count, function(n) {
      sub <- igraph::delete.vertices(graph, n)
      w - sum(1/rsfcNet::current_centr(sub))
    } )
    return(nodes)
  } else if (method=="regular") {
    w <- sum(1/igraph::closeness(graph, normalized = FALSE, weights = abs(E(graph)$weight)))
    count = igraph::vcount(graph)
    nodes <- sapply(1:count, function(n) {
      sub <- igraph::delete.vertices(graph, n)
      w - sum(1/igraph::closeness(sub, normalized = FALSE, weights = abs(E(sub)$weight)))
    } )
    return(nodes)
  }
}


#' Find the closeness vitality for a multiple graphs.
#'
#' Closeness vitality of a node is the change in the sum of distances between all node pairs when that node is removed from the network.
#' @param graphs A list of igraph objects or connectivity matrices.
#' @param method One of "regular" (the default) or "current"
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @param parallel Should multiple cores be used? Defaults to FALSE. If TRUE, progress bar is not displayed. This is normal.
#' @param cores How many cores should be used? Defaults to recommended 1 less than number of CPU cores.
#' @return A matrix containing the vitality scores for each node in each subject.
#' @author Brandon Vaughan
#' @export
#'
#' @details
#' Closeness vitality of a node is the change in the sum of distances between all node pairs when that node is removed from the network.
#' Either traditional closeness centrality can be used to define the sum of distances (Wiener Index) or the current closeness
#' centrality metric can be substituted with method="current". This function uses the absolute value of the weights, so if you require
#' only positive edges to be considered you must threshold them out of the network.
#'
#' The formula is fairly straightforward:
#'
#' First compute the Wiener Index for the entire graph, which is the sum of the reciprocals of the closeness
#' centrality for each node:
#'
#'
#' \eqn{W=\sum_{i=n}^N \frac{1}{ClC_n}}
#'
#'
#' Next, for each node, calculate the Wiener index again for a subgraph including all nodes but node_i. The vitality
#' of the node is simply the difference between the total Wiener index and the Wiener index for the subgraph.
#'
#'
#' \eqn{Vitality = W_{total} - W_{subgraph}}
#'
#'
#' The vitality of a node measures the increase in transport cost for the whole network
#' when a node is dropped. Negative values indicate a node can be excluded and decrease costs.
#' As the name suggests, the metric measures who in a network is "vital to the operation" and
#' who is an expense to the network (Brandes & Erlebach, 2005). Also see
#' \code{\link[rsfcNet:laplace_centr_mult]{laplacian centrality}}
#' for another measure of node centrality based on the impact to the network after deletion.
#'
#' @references
#' Brandes, U. & Erlebach, T. (2005). Network Analysis: Methodological Foundations, U.S. Government Printing Office.
#'
#' Brandes U., Fleischer D. (2005) Centrality Measures Based on Current Flow. In: Diekert V., Durand B. (eds) STACS 2005. STACS 2005. Lecture Notes in Computer Science, vol 3404. Springer, Berlin, Heidelberg
#'
#' @seealso
#' \code{\link[rsfcNet]{closeness_centr}}
#' \code{\link[rsfcNet]{current_centr}}
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#'
#' @examples
#'
#' **## Not run:**
#' vitality_mult(graphs)
#' ## End(**Not run**)
#'
#'
#'
vitality_mult = function(graphs, method="regular", col.names=NULL, row.names=NULL, parallel=TRUE, cores=NA){
if (parallel==FALSE){
  vit = pbapply::pbsapply(graphs, function(x) vitality(x, method=method))
  vit= t(vit)
  colnames(vit) = col.names
  rownames(vit) = row.names
  return(vit)
} else if (parallel==TRUE) {
  if (is.na(cores)==TRUE) {
    cl = as.integer(parallel::detectCores()-1)
    cl = parallel::makeCluster(cl)
  } else if (is.na(cores)==FALSE) {
    cl = as.integer(cores)
    cl = parallel::makeCluster(cl)
  }

  parallel::clusterEvalQ(cl, {
    library(parallel)
    library(igraph)
  })

  parallel::clusterExport(cl, "method", "graphs", envir=environment())

  vit = parallel::parLapply(cl, graphs, function(x) rsfcNet::vitality(x))
  parallel::stopCluster(cl)
  vit= t(vit)
  colnames(vit) = col.names
  rownames(vit) = row.names
  return(vit)
  }
}
