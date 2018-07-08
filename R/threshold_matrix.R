#' Threshold a correlation matrix
#'
#' This function is used to threshold a correlation matrix in order to prune weak connections.
#' @param c The correlation matrix.
#' @param q The quantile to threshold. Defaults to .95, such that all below the 95th quantile are set to zero.
#' @return The pruned correlation matrix.
#' @export
#' @author Brandon Vaughan
#' @examples
#' threshold_matrix(matrix, q=.90) # This should threshold such that only positive values are returned.
#' threshold_matrix(abs(matrix), q=.90) # Taking the absolute value of the matrix will allow negative correlations sufficiently strong to be included as an edge. This does not preserve signs.
#' threshold_matrix(abs(matrix), q=.90)*sign(matrix) This will restore the signs after thresholding if desired.
threshold_matrix = function(c, q=.95){
  threshold = quantile(c, q)
  c[which(c<threshold)] <- 0; c
}


