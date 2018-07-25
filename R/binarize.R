#' Binarize a correlation matrix.
#'
#' Binarize a correlation matrix.
#' @param c The correlation matrix.
#' @return A binary adjacency matrix
#' @export
#' @author Brandon Vaughan
#' @examples
#' binarize(graph)
#'
binarize = function (c, q = 0.9, method="quantile") {
  abs(sign(c))
}
