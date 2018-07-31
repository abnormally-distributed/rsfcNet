#' Calculate betweenness centrality for a list of graphs.
#'
#' This function is a convenience wrapper for igraph's betweenness centrality function and takes as an input a list of igraph objects.
#' @param graphs A list of igraph objects.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @param parallel Should multiple cores be used? Defaults to FALSE. If TRUE, progress bar is not displayed. This is normal.
#' @param cores How many cores should be used? Defaults to recommended 1 less than number of CPU cores.
#' @return A matrix of the betweenness centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' Betweenness centrality gives the number of times a given node lies in the shortest path between two other nodes.
#' Betweenness centrality may be less natural to interpret in rsfc networks (see Power et al 2013)
#' but remains a popular centrality metric. See for example Wang, Zuo, & He (2010).
#' It is calculated with the formula below:
#'
#' \eqn{b_i = \sum_{s \neq v \neq t}\frac{\sigma_{st}(i)}{\sigma_{st}}}
#'
#' @seealso
#' \code{\link[rsfcNet]{btwn_centr}}
#'
#' @examples
#' betweenness = btwn_centr_mult(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @references
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Power, J. D., Schlaggar, B. L., Lessov-Schlaggar, C. N., & Petersen, S. E. (2013). Evidence for hubs in human functional brain networks. Neuron, 79(4), 10.1016/j.neuron.2013.07.035. http://doi.org/10.1016/j.neuron.2013.07.035
#'
#' Wang, J., Zuo, X., & He, Y. (2010). Graph-Based Network Analysis of Resting-State Functional MRI. Frontiers in Systems Neuroscience, 4, 16. http://doi.org/10.3389/fnsys.2010.00016
#'
btwn_centr_mult = function(graphs,col.names=NULL, row.names=NULL, parallel=FALSE, cores=NA) {
  if (parallel==FALSE) {
  betweenness = pbapply::pblapply(graphs, function(x) igraph::betweenness(x, weights = abs(E(x)$weight), normalize=FALSE))
  betweenness = as.data.frame(t(data.frame(betweenness)))

  colnames(betweenness) = col.names
  rownames(betweenness) = row.names
  betweenness = as.matrix(betweenness)
  return(betweenness)
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

    betweenness = parallel::parLapply(cl, graphs, function(x) igraph::betweenness(x, weights = abs(E(x)$weight), normalize=FALSE))
    stopCluster(cl)
    betweenness = t(data.frame(betweenness))
    colnames(betweenness) = col.names
    rownames(betweenness) = row.names
    betweenness = as.matrix(betweenness)
    return(betweenness)

  }
}


#' Calculate betweenness centrality for a single graph.
#'
#' This function is a convenience wrapper for igraph's betweenness centrality function for consistency with rsfcNet's naming conventions.
#' @param graph A network as an igraph object or matrix.
#' @param weights a vector of edge weights. Defaults to the absolute value of the input network.
#' @return A vector of centrality scores for each node.
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' Betweenness centrality gives the number of times a given node lies in the shortest path between two other nodes.
#' Betweenness centrality may be less natural to interpret in rsfc networks (see Power et al 2013)
#' but remains a popular centrality metric. See for example Wang, Zuo, & He (2010).
#' It is calculated with the formula below:
#'
#' \eqn{b_i = \sum_{s \neq v \neq t}\frac{\sigma_{st}(i)}{\sigma_{st}}}
#'
#' @seealso
#' \code{\link[rsfcNet]{btwn_centr_mult}}
#'
#' @examples
#' betweenness = btwn_centr(graph)
#'
#' @references
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Power, J. D., Schlaggar, B. L., Lessov-Schlaggar, C. N., & Petersen, S. E. (2013). Evidence for hubs in human functional brain networks. Neuron, 79(4), 10.1016/j.neuron.2013.07.035. http://doi.org/10.1016/j.neuron.2013.07.035
#'
#' Wang, J., Zuo, X., & He, Y. (2010). Graph-Based Network Analysis of Resting-State Functional MRI. Frontiers in Systems Neuroscience, 4, 16. http://doi.org/10.3389/fnsys.2010.00016
#'
#'
btwn_centr = function(graph, weights = abs(E(graph)$weight)) {
  if (igraph::is.igraph(graph)==FALSE){
    graph = igraph::graph.adjacency(graph, mode="undirected", weighted=T, diag=FALSE)
  }
  igraph::betweenness(graph, weights = weights, normalize=FALSE)
}

