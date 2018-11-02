#' Obtain correlation matrices for each subject.
#'
#' Get regularized estimates of the correlation matrices using a variety of algorithms.
#' @param scrubbed The list of scrubbed time series matrices.
#' @param n The number of subjects in the list. Checks the global environment for "n" but may be assigned directly.
#' @param n.nodes The number of nodes. Checks the global environment for "n.nodes" but may be assigned directly.
#' @param method  One of "pearson", "partial", or "covar".
#' @return The connectivity matrices.
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' For methods "pearson", "partial", "covar", the estimation is done using the James-Stein based
#' methods described in Schaefer and Strimmer (2005) and the optimal shrinkage level calculated
#' for each subject using the method of Opgen-Rhein and Strimmer (2007). Respectively these give the
#' regularized pearson correlation matrix, shrinkage estimated partial correlation matrix, and regularized
#' covariance matrix, This depends on the corpcor package.The covariance option shouldn't be used directly
#' for network analysis, but can be used with whatever you need it for. As the number of observations increases
#' relative to the number of variables these all converge onto the empirical pearson, partial, and covariance
#' matrices respectively.
#'
#' @examples
#' \donttest{
#' cormats = get_cor_matrices(scrubbed_ts_list)
#' }
#'
#' @seealso
#' \code{\link[rsfcNet]{threshold_matrix}}
#' \code{\link[rsfcNet]{binarize}}
#'
#' @references
#'
#'
#' Opgen-Rhein, R., and K. Strimmer (2007). Accurate ranking of differentially expressed genes by a distribution-free shrinkage approach. Statist. Appl. Genet. Mol. Biol. 6:9. doi:10.2202/1544-6115.1252
#'
#' Schaefer, J., and K. Strimmer (2005). A shrinkage approach to large-scale covariance estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32. doi:10.2202/1544-6115.1175
#'
get_cor_matrices = function(scrubbed, method="partial", n=NULL, n.nodes=NULL) {

  if (method=="covar") {
    print("Estimating Regularized Covariance Matrices with Shrinkage.")
    covmats = pbapply::pblapply(scrubbed,corpcor::cov.shrink, verbose=FALSE)
    covmats = lapply(covmats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    return(covmats)
  } else if (method=="pearson") {
    print("Estimating Regularized Pearson Connectivity Matrices with Shrinkage.")
    cormats= pbapply::pblapply(scrubbed,corpcor::cor.shrink, verbose=FALSE)
    cormats= lapply(1:n, function(i) {diag(cormats[[i]]) <- 0.00000000 ; cormats[[i]]})
    cormats = lapply(cormats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    return(cormats)
  } else if (method=="partial") {
    print("Estimating Partial Connectivity Matrices With Shrinkage.")
    cormats= pbapply::pblapply(scrubbed,corpcor::pcor.shrink, verbose=FALSE)
    cormats= lapply(cormats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    cormats= lapply(1:n, function(i) {diag(cormats[[i]]) <- 0.00000000 ; cormats[[i]]})
    return(cormats)
   }
}




