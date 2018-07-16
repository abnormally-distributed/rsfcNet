#' Read Bad Time Points from Confound Matrix
#'
#' This identifies time points that must be scrubbed away. The function expects a confound matrix to be in the form provided by the fsl_motion_outliers function in the FSL software library.
#' @param confound_matrix The matrix of confounds.
#' @return A vector indicating which time points are in need of scrubbing.
#' @export
#' @examples
#' bad_ts(subj001_confoundm)
#' @author Brandon Vaughan
#' @references
#' https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLMotionOutliers
#'
bad_ts = function(confound_matrix) {which(rowSums(confound_matrix) > 0)}
