#' Obtain correlation matrices for each subject.
#'
#' This function is a convenience wrapper for various functions in the corpcor package. The default is the partial correlation matrix (method="partial") but the covariance matrix (method="covariance"), partial covariance (method="partial.var"), and regular Pearson correlation matrix (method="pearson") are also available options.  The estimation is done using the method described in Schafer and Strimmer (2005) and the optimal shrinkage level calculated for each subject using the method of Opgen-Rhein and Strimmer (2007).
#' @param scrubbed_ts_bulk The list of scrubbed time series matrices.
#' @param n The number of subjects in the list. Checks the global environment for "n" but may be assigned directly.
#' @param n.nodes The number of nodes. Checks the global environment for "n.nodes" but may be assigned directly.
#' @return The scrubbed time series matrices.
#' @export
#' @author Brandon Vaughan
#'
#' @examples
#' cormats = get_cor_matrices(scrubbed_ts_list)
#'
#' @references
#' Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of differentially expressed genes by a distribution-free shrinkage approach. Statist. Appl. Genet. Mol. Biol. 6:9. doi:10.2202/1544-6115.1252
#'
#'
#' Schaefer, J., and K. Strimmer. 2005. A shrinkage approach to large-scale covariance estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32. doi:10.2202/1544-6115.1175
#'
get_cor_matrices = function(scrubbed_ts_bulk, method="partial", n=NULL, n.nodes=NULL) {
  print("Estimating Connectivity Matrix with Shrinkage.")
  if (method=="covariance") {
    covmats = pbapply::pblapply(scrubbed_ts_bulk,corpcor::cov.shrink, verbose=FALSE)
    covmats= lapply(1:n, function(i) {diag(covmats[[i]]) <- 0; covmats[[i]]})
    covmats = lapply(covmats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    return(covmats)
  } else if (method=="pearson") {
    cormats= pbapply::pblapply(scrubbed_ts_bulk,corpcor::cor.shrink, verbose=FALSE)
    cormats= lapply(1:n, function(i) {diag(cormats[[i]]) <- 0; cormats[[i]]})
    cormats = lapply(cormats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    return(cormats)
  } else if (method=="partial") {
    cormats= pbapply::pblapply(scrubbed_ts_bulk,corpcor::pcor.shrink, verbose=FALSE)
    cormats= lapply(1:n, function(i) {diag(cormats[[i]]) <- 0; cormats[[i]]})
    cormats= lapply(cormats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    return(cormats)
  } else if (method=="partial.var"){
    covmats= pbapply::pblapply(scrubbed_ts_bulk,corpcor::pvar.shrink, verbose=FALSE)
    covmats= lapply(1:n, function(i) {diag(covmats[[i]]) <- 0; covmats[[i]]})
    covmats= lapply(covmats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    return(covmats)
  }
}
