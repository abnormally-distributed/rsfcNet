<<<<<<< HEAD
#' Obtain correlation matrices for each subject.
#'
#' Get regularized estimates of the correlation matrices using a variety of algorithms.
#' @param scrubbed The list of scrubbed time series matrices.
#' @param n The number of subjects in the list. Checks the global environment for "n" but may be assigned directly.
#' @param n.nodes The number of nodes. Checks the global environment for "n.nodes" but may be assigned directly.
#' @param method  One of "pearson", "partial", "covar", "parvar", "ridge", "elastic", or "adaglasso". See details.
#' @param cores How many cores for ridge or adaglasso
#' @param k how many cross validation folds to use in "adaglasso", "ridge", and "elastic". Defaults to 5.
#' @param alpha the penalty for method="elastic". 0 gives the ridge penalty (L2) and 1 gives the lasso penalty (L1). Defaults to .75.
#' @param pre.thresh If pre.thresh = TRUE (the default) correlations with an absolute value below .03 are set to zero.
#' @return The connectivity matrices.
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' For methods "pearson", "partial", "covar", and "parvar" the estimation is done using the James-Stein based
#' methods described in Schaefer and Strimmer (2005) and the optimal shrinkage level calculated
#' for each subject using the method of Opgen-Rhein and Strimmer (2007). Respectively these give the
#' regularized pearson correlation matrix, shrinkage estiamted partial correlation matrix, regularized
#' covariance matrix, and sihrinkage estimated partial variance matrix. This depends on the
#' corpcor package.The covariance and partial variance options shouldn't be used directly
#' for network analysis, but can be used with whatever you need them for. The option "ridge"
#' uses ridge regression (L2 penalty) with a cross-validated  penalty parameter to estimate
#' the partial correlation matrix. This method was developed by Kraemer, Schaefer, & Boulesteix (2009). This method
#' is much slower than the James-Stein type penalty used in Schaefer and Strimmer's algorithm,
#' but should yield more sparse results. For "adaglasso" the adaptive GLASSO of Zou (2006) method
#' recommended by Zhu & Cribben (2018) is used where an initial GLASSO  is used as a starting point
#' for the adaptive GLASSO procedure. For this method a list of lists is returned, where each
#' list contains both the best GLASSO solution and the best adaptive GLASSO
#' solution. Inspect both to see which is best suited to your needs. The ridge and adaptive GLASSO
#' methods depend on the pcor package. A final option provided by
#' this package is the elastic net method. This relies on the CorShrink package and allows the user to enter
#' an elastic net penalty to yield a compromise between a graphical ridge and graphical lasso estimate of the
#' partial correlation matrix.
#'
#'
#' Which method you should use depends on your goals. My recommendation is to start with
#' method=="partial" or method=="pearson" and leave the pre.thresh option on (the default).
#' Most likely this will produce a matrix sparse enough to work well with most functions.
#' If this is not sparse enough, consider thresholding further. Alternatively, use one of
#' the other methods. The ridge and adaptive GLASSO methods can be used to achieve greater sparsity.
#' The adaptive GLASSO method on the will return the GLASSO and adaptive GLASSO estimated matrices which are guaranteed to be sparse.
#' If the results are TOO sparse, however, the ridge or elastic methods should be attempted. The elastic net option
#' can be tried when the ridge isn't sparse enough and the adaptive lasso or lasso outputs are too sparse.
#' The elastic method also is a little faster than the ridge and adaptive GLASSO options. Highly sparse networks
#' such as that obtained with GLASSO may be desirable if you want to create a binary matrix by setting nonzero
#' edges to 1 later, but if you want to work with weighted edges the ridge and James-Stein estimators are
#' probably more ideal.
#'
#' For methods "pearson", "partial", "ridge", and "elastic" the option pre.thresh sets the smallest
#' correlations (those with an absolute value below 0.015, or 0.1 for method="pearson")  to zero. Correlations this small are likely
#' just noise, and are of no interest to estimating a graphical model. These set default values were chosen
#' based on candidate values when testing different pre-thresholds on private data. These values seem to work
#' well to prevent disconnected nodes (which breaks certain graph metrics) while also inducing enough sparsity for functions that work poorly
#' on full-edged dense networks to be able to calculate results. Further thresholding can be applied with the threshold_matrix function. This feature can be turned
#' off simply by setting pre.thresh=FALSE.
#'
#' @examples
#' cormats = get_cor_matrices(scrubbed_ts_list)
#'
#' @seealso
#' \code{\link[rsfcNet]{threshold_matrix}}
#' \code{\link[rsfcNet]{threshold_and_binarize}}
#'
#' @references
#'
#' Kraemer, N., Schaefer, J., and Boulesteix, A. (2009) "Regularized Estimation of Large-Scale Gene Regulatory Networks using Gaussian Graphical Models", BMC Bioinformatics, 10:384
#'
#' Opgen-Rhein, R., and K. Strimmer (2007). Accurate ranking of differentially expressed genes by a distribution-free shrinkage approach. Statist. Appl. Genet. Mol. Biol. 6:9. doi:10.2202/1544-6115.1252
#'
#' Schaefer, J., and K. Strimmer (2005). A shrinkage approach to large-scale covariance estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32. doi:10.2202/1544-6115.1175
#'
#' Stephens, M. (2017). False discovery rates: a new deal, Biostatistics, Volume 18, Issue 2, Pages 275–294, https://doi.org/10.1093/biostatistics/kxw041
#'
#' Zou, H. (2006) "The Adaptive Lasso and its Oracle Property", Journal of the American Statistical Association. 101 (476): 1418-1429.
#'
#' Zhu, Y., & Cribben, I. (2018). Sparse Graphical Models for Functional Connectivity Networks: Best Methods and the Autocorrelation Issue. Brain Connectivity, 8(3), 139-165. doi:10.1089/brain.2017.0311
#'
get_cor_matrices = function(scrubbed, method="partial", n=NULL, n.nodes=NULL, cores=NA, k=5, alpha=.75, pre.thresh=TRUE) {

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
    if (pre.thresh==TRUE){
     cormats = lapply(1:n, function(i) {cormats[[i]][which(abs(cormats[[i]]) < 0.1)] <- 0.00000000 ; cormats[[i]]})
    }
    return(cormats)
  } else if (method=="partial") {
    print("Estimating Partial Connectivity Matrices With Shrinkage.")
    cormats= pbapply::pblapply(scrubbed,corpcor::pcor.shrink, verbose=FALSE)
    cormats= lapply(cormats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    cormats= lapply(1:n, function(i) {diag(cormats[[i]]) <- 0.00000000 ; cormats[[i]]})
    if (pre.thresh==TRUE){
    cormats = lapply(1:n, function(i) {cormats[[i]][which(abs(cormats[[i]]) < 0.015)] <- 0.00000000 ; cormats[[i]]})
    }
    return(cormats)
  } else if (method=="parvar"){
    print("Estimating Partial Variance Matrices With Shrinkage.")
    covmats= pbapply::pblapply(scrubbed,corpcor::pvar.shrink, verbose=FALSE)
    covmats= lapply(covmats,function(i) matrix(i, nrow=n.nodes, ncol=n.nodes))
    return(covmats)
  } else if (method=="ridge") {
    cores = parallel::detectCores()
    if (cores==1) {
      print("Ridge estimation requires a multicore processor")
    } else {
      cl = as.integer(parallel::detectCores()-1)
      cl = parallel::makeCluster(cl)
      parallel::clusterExport(cl, "k", envir=environment())
      parallel::clusterEvalQ(cl, {
        library(parcor)
        library(parallel)
      })
      print("Estimating Partial Connectivity Matrices With Cross Validated Ridge Penalty.")
      print("This may take a long time. For a faster method choose another method.")
      ridge_parcor = parallel::parLapply(cl,scrubbed, function(x) ridge.net(x, k=k, verbose=FALSE)$pcor)
      parallel::stopCluster(cl)
      cormats= lapply(1:n, function(i) {diag(ridge_parcor[[i]]) <- 0.00000000 ;  ridge_parcor[[i]]})
      if (pre.thresh==TRUE){
        cormats = lapply(1:n, function(i) {cormats[[i]][which(abs(cormats[[i]]) < 0.015)] <- 0.00000000 ; cormats[[i]]})
      }
      return(cormats)
      }
      } else if (method=="adaglasso") {
      cores = parallel::detectCores()
      if (cores==1) {
        print("Adaptive GLASSO estimation requires a multicore processor")
      } else {
        cl = as.integer(parallel::detectCores()-1)
        cl = parallel::makeCluster(cl)
        parallel::clusterExport(cl, "k", envir=environment())
        parallel::clusterEvalQ(cl, {
          library(parcor)
          library(parallel)
        })
        print("Estimating Partial Connectivity Matrices With Cross Validated Adaptive GLASSO")
        print("This will take a very long time. For a faster method choose another method.")
        print("Also consider leaving this to run overnight or running on a server.")
        ada_parcor = parallel::parLapply(cl,scrubbed, function(x) adalasso.net(x, k=k, verbose=FALSE, intercept=FALSE))
        parallel::stopCluster(cl)
        cormats= lapply(1:n, function(i) {diag(ada_parcor[[i]]) <- 0.00000000 ;  ada_parcor[[i]]})
        return(cormats)
      }
      } else if (method=="elastic") {
        cores = parallel::detectCores()
        if (cores==1) {
          print("Elastic Net estimation requires a multicore processor")
        } else {
          cl = as.integer(parallel::detectCores()-1)
          cl = parallel::makeCluster(cl)
          parallel::clusterExport(cl, c("k", "alpha"), envir=environment())
          parallel::clusterEvalQ(cl, {
            library(CorShrink)
            library(parallel)
          })
          print("Estimating Partial Connectivity Matrices With Cross Validated Elastic Net")
          print("This could take a while. For a faster method choose another method.")
          ada_parcor = parallel::parLapply(cl,scrubbed, function(x) CorShrink::pCorShrinkData(x, glmnet_alpha =  alpha, glmnet_nfolds=k))
          parallel::stopCluster(cl)
          cormats= lapply(1:n, function(i) {diag(ada_parcor[[i]]) <- 0.00000000 ;  ada_parcor[[i]]})
          if (pre.thresh==TRUE){
            cormats = lapply(1:n, function(i) {cormats[[i]][which(abs(cormats[[i]]) < 0.015)] <- 0.00000000 ; cormats[[i]]})
          }
          return(cormats)
        }
      }
}




=======
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
>>>>>>> 48e30e2d1843035d6169a05c358a6b61b43f39fb
