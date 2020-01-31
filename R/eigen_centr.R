
#' Calculate eigenvector centrality for a single graph.
#'
#' @description This function calculates the eigenvector centrality for a single graph.
#'
#' @param graph A network as an igraph object or matrix.
#' @param normalize how the normalization of eigenvector centrality should be treated. The options are as follows: \cr
#' \cr
#' If 0, no normalization is applied.  \cr
#' If 1, the absolute value of the principal eigenvector is used as in the MATLAB Brain Connectivity Toolbox (Rubinov & Sporns, 2010),
#' but no normalization is applied.  \cr
#' If 2, the principal eigenvector is taken as an absolute value, and the result is divided
#' by the maximum value such that the largest eigenvector centrality is 1.  \cr
#' If 3, the ranks of the raw eigenvector are taken, and the ranks divided by the maximum rank such that the largest eigenvector centrality is 1. \cr
#' If 4, eigenvector entries less than zero are set to zero, and then the scores are divided by the maximum value such that the largest eigenvector centrality is 1. \cr
#' If 5, the raw ranks of the principal eigenvector entries are returned. \cr
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
eigen_centr = function(graph, normalize = c(0, 1, 2, 3, 4, 5)) {
  if (is.igraph(graph)=="TRUE"){
    graph = as.matrix(igraph::as_adjacency_matrix(graph, edges = FALSE, attr = "weight", sparse=TRUE)) }

  diag(graph) <- 1
  EVC <-eigen(graph) # get eigenvalues and eigenvectors
  eigen_centrality <- EVC$vectors[,1] # get the eigenvector associated with the largest eigenvalue
  normalize <- as.character(normalize)
  normalize <- match.arg(normalize, as.character(c(0,1,2,3,4,5)))
  normalize <- as.numeric(noquote(normalize))
  if (normalize == 0){ ## none
    if (all(eigen_centrality < 0)) eigen_centrality <- eigen_centrality * -1
    eigen_centrality <- eigen_centrality
  }
  else if (normalize == 1){ ## absval
    eigen_centrality <- abs(eigen_centrality)
  }
  else if (normalize == 2){ ## max_absval
    eigen_centrality <- abs(eigen_centrality)/max(abs(eigen_centrality))
  }
  else if (normalize == 3){ ## max
    if (all(eigen_centrality < 0)) eigen_centrality <- eigen_centrality * -1
    eigen_centrality <- rank(eigen_centrality) / max(rank(eigen_centrality))
  }
  else if (normalize == 4){ ## thresh
    if (all(eigen_centrality < 0)) eigen_centrality <- eigen_centrality * -1
    eigen_centrality[which(eigen_centrality < 0)] <- 0; eigen_centrality <- eigen_centrality/max(eigen_centrality)
  }
  else if (normalize == 5){ ## rank
    if (all(eigen_centrality < 0)) eigen_centrality <- eigen_centrality * -1
    eigen_centrality <- rank(eigen_centrality)
  }
  return(eigen_centrality)
}


#' Calculate eigenvector centrality for multiple graphs.
#'
#' @description This function calculates the eigenvector centrality for multiple graphs.
#'
#' @param graphs a list of igraph objects.
#' @param normalize how the normalization of eigenvector centrality should be treated. The options are as follows: \cr
#' \cr
#' If 0, no normalization is applied.  \cr
#' If 1, the absolute value of the principal eigenvector is used as in the MATLAB Brain Connectivity Toolbox (Rubinov & Sporns, 2010),
#' but no normalization is applied.  \cr
#' If 2, the principal eigenvector is taken as an absolute value, and the result is divided
#' by the maximum value such that the largest eigenvector centrality is 1.  \cr
#' If 3, the ranks of the raw eigenvector are taken, and the ranks divided by the maximum rank such that the largest eigenvector centrality is 1. \cr
#' If 4, eigenvector entries less than zero are set to zero, and then the scores are divided by the maximum value such that the largest eigenvector centrality is 1. \cr
#' If 5, the raw ranks of the principal eigenvector entries are returned. \cr
#' @return A vector of centrality scores for each node.
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
eigen_centr_mult = function(graphs,normalize = c(0, 1, 2, 3, 4, 5), col.names = NULL, row.names = NULL) {
  eigen.centrality = pbapply::pbsapply(graphs, function(x) rsfcNet::eigen_centr(x, normalize))
  eigen.centrality = t(eigen.centrality)
  colnames(eigen.centrality) = col.names
  rownames(eigen.centrality) = row.names
  eigen.centrality = as.matrix(eigen.centrality)
  return(eigen.centrality)
}
