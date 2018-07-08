#' Imports Time Series or Confound Files from the list created with get_file_list.
#'
#' This function loads the files into the global environment as a list of matrices. The function expects each file to be a .csv.
#' @param file.list The list of files.
#' @return The time series matrix or confound matrix for each subject.
#' @export
#' @author Brandon Vaughan
#'
#' @examples
#' ts = import_from_list(ts_file_locations)
#'
#' confounds = import_from_list(confound_file_locations)
#'
import_from_list = function(file.list) {do.call(list,lapply(file.list, function(x) read.csv(file=x,sep= "",header = TRUE)))}
