#' Calculate eigenvector centrality for all nodes for multiple subjects.
#'
#' This function is a convenience wrapper for igraph's eigenvector centrality function and takes as an input a list of igraph objects. ECV is calculated on the absolute values of the correlation matrix per Lohmann et al, 2010. ECV is standardized subject-wise to 0-1 range for consistency.
#' @param graphs a list of igraph objects.
#' @param col.names The names of each column (node labels). Checks global environment for "colnames" but may be assigned directly.
#' @param row.names The names of each row (subject). Checks global environment for "rownames" but may be assigned directly.
#' @return A matrix of the eigenvector centralities of each node for each subject.
#' @export
#' @author Brandon Vaughan
#' @examples
#' eigencent = eigen_centr(graphs,row.names = subj_numbers,col.names = node_labels)
#'
#' @references
#' Lohmann, G., Margulies, D. S., Horstmann, A., Pleger, B., Lepsien, J., Goldhahn, D., … Turner, R. (2010). Eigenvector Centrality Mapping for Analyzing Connectivity Patterns in fMRI Data of the Human Brain. PLoS ONE, 5(4), e10232. doi:10.1371/journal.pone.0010232
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4

eigen_centr = function(graphs, col.names=.GlobalEnv$colnames, row.names=.GlobalEnv$rownames) {
  eigen.centrality = pbapply::pbsapply(graphs, function(x) igraph::eigen_centrality(x,directed=FALSE,scale = TRUE,weights=abs(E(x)$weight))$vector)
  eigen.centrality = t(eigen.centrality)

  colnames(eigen.centrality) = col.names
  rownames(eigen.centrality) = row.names
  eigen.centrality = as.matrix(eigen.centrality)
  return(eigen.centrality)
}
