#' Calculate betweenness centrality for all nodes for multiple subjects.
#'
#' This function is a convenience wrapper for igraph's betweenness centrality function and takes as an input a list of igraph objects. Betweenness centrality gives the number of times a given node lies in the shortest path between two other nodes. Betweenness centrality may be less natural to interpret in rsfc networks (see Power et al 2013) but remains a popular centrality metric. See for example Wang, Zuo, & He (2010).
#' @param graphs a list of igraph objects.
#' @param col.names The names of each column (node labels). Checks global environment for "colnames" but may be assigned directly.
#' @param row.names The names of each row (subject). Checks global environment for "rownames" but may be assigned directly.
#' @return A matrix of the betweenness centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#' @examples
#' betweenness = btwn_centr(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @references
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Power, J. D., Schlaggar, B. L., Lessov-Schlaggar, C. N., & Petersen, S. E. (2013). Evidence for hubs in human functional brain networks. Neuron, 79(4), 10.1016/j.neuron.2013.07.035. http://doi.org/10.1016/j.neuron.2013.07.035
#'
#' Wang, J., Zuo, X., & He, Y. (2010). Graph-Based Network Analysis of Resting-State Functional MRI. Frontiers in Systems Neuroscience, 4, 16. http://doi.org/10.3389/fnsys.2010.00016


btwn_centr = function(graphs,col.names=.GlobalEnv$colnames, row.names=.GlobalEnv$rownames) {
  betweenness = pbapply::pblapply(graphs, function(x) igraph::betweenness(x, weights = abs(E(x)$weight), normalize=FALSE))
  betweenness = as.data.frame(t(data.frame(betweenness)))

  colnames(betweenness) = col.names
  rownames(betweenness) = row.names
  betweenness = as.matrix(betweenness)
  return(betweenness)
}
