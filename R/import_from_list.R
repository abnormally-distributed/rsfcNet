#' Imports Time Series or Confound Files from the list created with get_file_list.
#'
#' This function loads the files into the global environment as a list of matrices. The function expects each file to be in a .csv or .txt format.
#' @param file.list The list of files.
#' @param sep The separator used to separate entries. Defaults to "" since this is the output format from fsl.
#' @param header Are the cells in the top row of the spreadsheet colum names? Must be TRUE or FALSE.
#' @param matrix Should the data be imported directly to a matrix? Defaults to TRUE.
#'
#' @return The time series matrix or confound matrix for each subject.
#' @export
#'
#'
#' @author Brandon Vaughan
#'
#' @examples
#' ts = import_from_list(ts_file_locations, header=TRUE)
#'
#' confounds = import_from_list(confound_file_locations, header=FALSE)
#'

import_from_list = function (file.list, header, sep = "", matrix = TRUE)
{

  read_matrix <-
    function(file, header = FALSE, sep = "", skip = 0)
    {
      row.lens <- count.fields(file, sep = sep, skip = skip)
      if(any(row.lens != row.lens[1]))
        stop("number of columns is not constant")
      if(header) {
        nrows <- length(row.lens) - 1
        ncols <- row.lens[2]
        col.names <- scan(file, what = "", sep = sep, nlines = 1,
                          quiet = TRUE, skip = skip)
        x <- scan(file, sep = sep, skip = skip + 1, quiet = TRUE)
      }
      else {
        nrows <- length(row.lens)
        ncols <- row.lens[1]
        x <- scan(file, sep = sep, skip = skip, quiet = TRUE)
        col.names <- NULL
      }
      x <- as.double(x)
      if(ncols > 1) {
        dim(x) <- c(ncols,nrows)
        x <- t(x)
        colnames(x) <- col.names
      }
      else if(ncols == 1)
        x <- as.vector(x)
      else
        stop("wrong number of columns")
      return(x)
    }



  if (matrix == TRUE) {
    list = lapply(file.list, function(x) read_matrix(file = x,
                                                     sep = sep, header = header))
    return(list)
  }
  else {
    list = lapply(file.list, function(x) read.csv(file = x,
                                                  sep = sep, header = header))
    return(list)
  }
}

