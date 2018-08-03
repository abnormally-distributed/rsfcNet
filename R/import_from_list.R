#' Imports Time Series or Confound Files from the list created with get_file_list.
#'
#' This function loads the files into the global environment as a list of matrices. The function expects each file to be in a .csv or .txt format.
#' @param file.list The list of files.
#' @param sep The separator used to separate entries. Defaults to "" since this is the output format from fsl.
#' @param header Are there column headers? Defaults to header = FALSE.
#' @param numeric Should all columnns be coercred to numeric if possible?
#'
#' @return The time series matrix or confound matrix for each subject.
#' @export
#'
#'
#' @author Brandon Vaughan
#'
#' @examples
#' ts = import_from_list(ts_file_locations)
#'
#' confounds = import_from_list(confound_file_locations)
#'
import_from_list = function(file.list, header = FALSE, sep="", numeric=TRUE){
  list = do.call(list,lapply(file.list, function(x) read.csv(file=x,sep= sep,header = header)))
  if (numeric==TRUE) {
    list = lapply(list, function(m) {apply(m, 2, function (c) as.numeric(c))})
    return(list)
  } else {
    list = lapply(list, function(m) {as.matrix(m)})
    return(list)
  }
}
