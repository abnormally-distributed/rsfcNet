#' Average together a list of connectivity matrices
#'
#' Average together a list of correlation matrices to analyze the typical network properties of a given group (ie, teens vs adults)
#' @param cormats A list of correlation matrices
#' @param n The sample size.
#' @param which can be "graph" to get an igraph object, "matrix" for a correlation matrix, or "both" for a list of both. Defaults to graph.
#' @return A matrix, igraph object, or list.
#' @export
#' @author Brandon Vaughan
#' @examples
#' \donttest{
#' cormats = get_cor_matrices(scrubbed_ts_list) #Get the correlation matrices first
#' graphs = get_average_graph(cormats)
#' }
#'
get_average_graph = function(cormats, n=NULL, which="graph"){
  n=length(cormats)
  average_matrix = lapply(1:n, function(i) {diag(cormats[[i]]) <- 0; cormats[[i]]})
  average_matrix = lapply(1:length(average_matrix), function(m) psych::fisherz(average_matrix[[m]]))
  average_matrix = apply(simplify2array(average_matrix),c(1,2), median)
  average_matrix = psych::fisherz2r(average_matrix)
  diag(average_matrix) <- 0

  average_graph = graph.adjacency(average_matrix, mode="undirected", weighted=T,diag = TRUE)
  if(which=="graph") {
    return(average_graph)
  } else if (which=="matrix") {
    return(average_matrix)
  } else if (which=="both"){
    return(list(average_matrix=average_matrix,average_graph=average_graph))
  }
}
