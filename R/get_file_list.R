#' Imports List of Time Series or Confound Files
#'
#' This function stores the file locations of each subject's file in a list.
#' @param path The location of the files.
#' @param full.names TRUE gets full path. FALSE retrieves only the names of the files. Defaults to TRUE.
#' @return The list of file locations or names.
#' @export
#' @author Brandon Vaughan
#' @examples
#' ts_file_locations = get_file_list(path="C:/Users/Yourname/Documents/brainstudy/ts/")
#'
#' confound_file_locations = get_file_list(path="C:/Users/Yourname/Documents/brainstudy/confoundm/")
#'
get_file_list = function(path,full.names="TRUE") {list.files(path=path,full.names = full.names)}
