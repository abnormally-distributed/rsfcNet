#' Calculate page_rank centrality for all nodes for multiple subjects.
#'
#' This function is a convenience wrapper for igraph's page_rank function.
#' @param graphs a list of igraph objects.
#' @param prior.probs = an optional list of vectors containing weights for the prior probability of each node being visited.
#' @param damp = damping factor determining the probability that a random walk will continue from a given node. Defaults to .85
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#'
#' @return A matrix of the page_rank centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#'
#' @details The 'pagerank' centrality type results from a random walk of the network. At each node in the graph,
#' the next node is jumped to with probability p from the neighbors of node i.
#' Otherwise, or when a node has no successors, the next node is chosen from all nodes.
#' The centrality score is the average time spent at each node during the random walk. Pagerank assumes
#' equal probability for each node, but prior probabilities can be supplied.
#'
#' @examples
#' page_rank = pagerank_centr(graphs)
#'
#' @seealso
#' \code{\link[rsfcNet]{eigen_centr}}
#'
#' @references
#'
#' Brin, S., Page, L. (1998). The Anatomy of a Large-Scale Hypertextual Web Search Engine. Proceedings of the 7th World-Wide Web Conference, Brisbane, Australia.
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
pagerank_centr = function (graphs, damp=.85, prior.probs="NULL", col.names = NULL, row.names = NULL)
{

  if (isTRUE(prior.probs=="NULL")==FALSE) {
    prior.probs =  lapply(prior.probs, function(x) {
      x <- as.matrix(x)
      minAttr=apply(x, 2, min)
      maxAttr=apply(x, 2, max)
      x <- sweep(x, 2, minAttr, FUN="-")
      x=sweep(x, 2,  maxAttr-minAttr, "/")
      return (x)})
    pagerank_centr = pbapply::pbsapply(1:as.numeric(length(graphs)), function(x) igraph::page_rank(graphs[[x]], directed=FALSE, damping=damp, personalized=prior.probs[[x]], weights= abs(E(graphs[[x]])$weight))$vector)

  } else {
    pagerank_centr = pbapply::pbsapply(graphs, function(x) igraph::page_rank(x, directed=FALSE, damping=damp, weights= abs(E(x)$weight))$vector)
  }

  pagerank_centr = t(pagerank_centr)
  colnames(pagerank_centr) = col.names
  rownames(pagerank_centr) = row.names
  pagerank_centr = as.matrix(pagerank_centr)
  return(pagerank_centr)
}
