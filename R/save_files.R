#' Save a list of files
#'
#' Save a list of matrices to the .csv or .txt file format.
#'
#' @param files a list of matrices with numeric data.
#' @param path a path to save files in. Required, but if left blank saves to current working directory.
#' @param names a vector of names for files. Required.
#' @param file.type one of ".csv", ".tsv" or ".txt". Defaults to .csv
#' @param sep the seperator to use for separating entries across columns. One of "," for comma separated values, " " for space separated values, or "'\\'t" for tab separated values. Defaults to ","
#' @param row.names Rowwise names to be used as the first column if desired. Defaults to none (row.names=FALSE).
#' @param col.names Column names to be used as the top row if desired. Defaults to none (col.names=FALSE).
#' @return returns nothing, but saves files to a directory.
#' @export
#' @author Brandon Vaughan
#' @examples
#' **## Not run:**
#' save_files(cormats, subj.names, path="C:/Users/YourName/Documents/R/")
#' ## End(**Not run**)
save_files =
  function (files, names, path = paste0(getwd(), "/") , file.type = ".csv",  sep = ",", row.names = FALSE, col.names = FALSE)
  {
    filenames <- paste0(names, file.type)
    files2 = lapply(files, function(x) as.matrix(x))
    for (i in 1:length(files)) {
      outname <- paste(path, filenames[i], sep = "")
      write.table(files2[[i]], outname, quote = F, row.names = row.names,
                  col.names = col.names, sep = sep)
    }
  }

