#' Used in scrub_time_series function.
#'
#' This function is used to scrub the ith matrix in a list. For use in scrub_time_series function.
#' @param i The time series matrix
#' @return The scrubbed time series matrix.
#' @export
#' @author Brandon Vaughan
#' @examples
#' See scrub_time_series.
scrub_matrix =  function(i){ts_file_bulk[[i]][-as.vector(badrows[[i]]),] }
