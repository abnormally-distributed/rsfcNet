#' Threshold a correlation matrix
#'
#' This function is used to threshold a correlation matrix in order to prune weak connections.
#' @param c The network either as connectivity matrix or as an igraph object.
#' @param q.thresh The quantile to threshold. Defaults to .90, such that all below the 90th quantile are set to zero.
#' @param thresh The number (taken to be an absolute value) below which to threshold (ie, correlation or cohen's d , etc). No default. Recommend .10 for correlation or .20 for cohen's d.
#' @return The pruned correlation matrix.
#' @export
#' @author Brandon Vaughan
#' @examples
#' threshold_matrix(matrix, q=.90) # This should threshold such that only positive values are returned.
#' threshold_matrix(abs(matrix), q=.90) # Taking the absolute value of the matrix will allow negative correlations sufficiently strong to be included as an edge. This does not preserve signs.
#' threshold_matrix(abs(matrix), q=.90)*sign(matrix) This will restore the signs after thresholding if desired.

threshold_matrix = function(c, method="quantile", q.thresh=.90, thresh){

  if (is.igraph(c)=="TRUE"){
  c = as.matrix(igraph::as_adjacency_matrix(c, edges = FALSE, attr = "weight", sparse=TRUE)) }
  if (method=="quantile") {
    threshold = quantile(abs(c), q.thresh)
    c[which(c<threshold)] <- 0; c
  }

}

#' Transform a correlation matrix into cohen's d or z-scores and back
#'
#' This function is used to transform correlation matrices into other units such as Cohen's d or Fisher's z.
#' @param c The network as a matrix
#' @param method one of "cohens.d", "undo.d"", "fishers.z" or "undo.d""
#' @return The transformed correlation matrix.
#' @export
#' @author Brandon Vaughan
#' @examples
#' transform_matrix(matrix, method="cohens.d") # This should threshold such that only positive values are returned.

transform_matrix = function(c, method="cohens.d") {
  if (method=="cohens.d") {
    2 * c/sqrt(1 - c^2)
  } else if (method=="undo.d"){
    c/sqrt(c^2 + 4)
  } else if (method=="fishers.z"){
      0.5 * log((1 + c)/(1 - c))
  } else if (method=="undo.z") {
    (exp(2 * c) - 1)/(1 + exp(2 * c))
  }
}



