#' Binarize a correlation matrix.
#'
#' Threshold a correlation matrix and convert it to an un-weighted (binary) adjacency matrix.
#' @param c The correlation matrix.
#' @param q The quantile to threshold. Defaults to .90, such that all below the 90th quantile are set to zero.
#' @return A binary adjacency matrix
#' @export
#' @author Brandon Vaughan
#' @examples
#' threshold_and_binarize(matrix, q=.90) # This should threshold such that only positive values are returned.
#' threshold_and_binarize(abs(matrix), q=.90) # Taking the absolute value of the matrix will allow negative correlations sufficiently strong to be included as an edge.
threshold_and_binarize = function(c, q){
  if (is.igraph(c)=="TRUE"){
  c = as.matrix(igraph::as_adjacency_matrix(c, edges = FALSE, attr = "weight", sparse=TRUE)) }
  m = threshold_matrix(c,q)
  abs(sign(m))
}
