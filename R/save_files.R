#' Save a list of numeric matrices
#'
#' Save a list of numeric matrices to the .csv file format.
#'
#' @param files a list of matrices with numeric data.
#' @param path a path to save files in. Required, but if left blank saves to current working directory.
#' @param names a vector of names for files. Required.
#' @return returns nothing, but saves files to a directory.
#' @export
#' @author Brandon Vaughan
#' @examples
#' **## Not run:**
#' save_files(cormats, subj.names, path="C:/Users/YourName/Documents/R/")
#' ## End(**Not run**)
save_files = function(files, names, path=getwd()) {
  filenames <- paste0(names, ".csv")
  files2 = lapply(files, function(x) as.matrix(x))
   for (i in 1:length(files)){
    outname <- paste(path, filenames[i], sep= "")
    write.table(files2[[i]], outname, quote= F, row.names = FALSE, col.names=FALSE, sep=",")
  }
}
