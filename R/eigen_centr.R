
#' Calculate eigenvector centrality for a single graph.
#'
#' This function calculates the eigenvector centrality for a single graph. This function
#' is essentially a port of the matlab script from the Brain Connectivity Toolbox (Rubinov & Sporns, 2010).
#'
#' @param graphs A network as an igraph object or matrix.
#' @return A vector of centrality scores for each node.
#' @export
#' @author Brandon Vaughan
#'
#' @details
#'
#' The eigenvector centrality is the eigenvector for each node corresponding to the largest eigenvalue of the matrix.
#'
#' Eigenvector centrality differs from degree or strength.
#' A node with many connections does not necessarily have a high eigenvector
#' centrality. For example, a node may have many very weak connections that yield a large value for strength/degree.
#' Likewise, a node with high eigenvector centrality may have few connections but be well connected to a small number
#' of important nodes.
#'
#' @examples
#'
#' \donttest{
#' eigencent = eigen_centr(graph)
#' }
#'
#' @seealso
#' \code{\link[rsfcNet]{leverage_centr}}
#' \code{\link[rsfcNet]{eigen_centr_mult}}
#'
#' @references
#'
#' \href{https://sites.google.com/site/bctnet/measures/list}{Brain Connectivity Toolbox}
#'
#' Lohmann, G., Margulies, D. S., Horstmann, A., Pleger, B., Lepsien, J., Goldhahn, D., … Turner, R. (2010). Eigenvector Centrality Mapping for Analyzing Connectivity Patterns in fMRI Data of the Human Brain. PLoS ONE, 5(4), e10232. doi:10.1371/journal.pone.0010232
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Rubinov, M. , Sporns, O. (2010) Complex network measures of brain connectivity: Uses and interpretations. NeuroImage 52:1059-69.
#'
eigen_centr = function(graph) {
  if (is.igraph(graph)=="TRUE"){
    graph = as.matrix(igraph::as_adjacency_matrix(graph, edges = FALSE, attr = "weight", sparse=TRUE)) }

  diag(graph) <- 0

  EVC <-eigen(graph) # get eigenvalues and eigenvectors
  eigen_centrality <- EVC$vectors[,1] # get the eigenvector associated with the largest eigenvalue
  eigen_centrality <- abs(eigen_centrality) #absolute value since the signs are arbitrary
  return(eigen_centrality)
}


#' Calculate eigenvector centrality for multiple graphs.
#'
#' This function calculates the eigenvector centrality for multiple graphs. This function
#' is essentially a port of the matlab script from the Brain Connectivity Toolbox (Rubinov & Sporns, 2010).
#'
#' @param graphs a list of igraph objects.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @return A matrix of the eigenvector centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' The eigenvector centrality is the eigenvector for each node corresponding to the largest eigenvalue of the matrix.
#'
#' Eigenvector centrality differs from degree or strength.
#' A node with many connections does not necessarily have a high eigenvector
#' centrality. For example, a node may have many very weak connections that yield a large value for strength/degree.
#' Likewise, a node with high eigenvector centrality may have few connections but be well connected to a small number
#' of important nodes.
#'
#' @examples
#'
#' \donttest{
#' eigencent = eigen_centr_mult(graphs,row.names = subj_numbers,col.names = node_labels)
#' }
#'
#' @seealso
#' \code{\link[rsfcNet]{eigen_centr}}
#' \code{\link[rsfcNet]{leverage_centr_mult}}
#'
#' @references
#'
#' \href{https://sites.google.com/site/bctnet/measures/list}{Brain Connectivity Toolbox}
#'
#' Lohmann, G., Margulies, D. S., Horstmann, A., Pleger, B., Lepsien, J., Goldhahn, D., … Turner, R. (2010). Eigenvector Centrality Mapping for Analyzing Connectivity Patterns in fMRI Data of the Human Brain. PLoS ONE, 5(4), e10232. doi:10.1371/journal.pone.0010232
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Rubinov, M. , Sporns, O. (2010) Complex network measures of brain connectivity: Uses and interpretations. NeuroImage 52:1059-69.
#'
eigen_centr_mult = function(graphs, col.names = NULL, row.names = NULL) {
  eigen.centrality = pbapply::pbsapply(graphs, function(x) rsfcNet::eigen_centr(x))
  eigen.centrality = t(eigen.centrality)
  colnames(eigen.centrality) = col.names
  rownames(eigen.centrality) = row.names
  eigen.centrality = as.matrix(eigen.centrality)
  return(eigen.centrality)
}
