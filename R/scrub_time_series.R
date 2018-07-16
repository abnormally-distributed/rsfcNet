#' "Scrubs" Bad Time Points from Time Series Matrix
#'
#' This function takes a list of time series files and a list of confound matrices and returns the scrubbed time series files. The default method is "censor". Censoring takes the preceding time point and the following two time points of each bad time point in the confound matrices and deletes them as described in Power et al, 2012. If there are relatively few bad time points in the subjects method="interpolate" may be tried.
#' @param ts_file_bulk list of time series matrices.
#' @param confound_file_bulk list of confound files for each subject
#' @param n The The number of subjects in the list.
#' @param n.nodes The The number of nodes in the network.
#' @param method method="censor" uses censoring/scrubbing while interpolation fills in missing time points column-wise with smooth polynomials.
#' @return The time series matrix for each subject.
#' @export
#' @author Brandon Vaughan
#' @examples
#' scrubbed_ts = scrub_time_series(ts, confounds)
#' @references
#' Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142–2154. http://doi.org/10.1016/j.neuroimage.2011.10.018
#'
#' Moritz S, Bartz-Beielstein T (2017). “imputeTS: Time Series Missing Value Imputation in R.” The R Journal, 9(1), 207–218. https://journal.r-project.org/archive/2017/RJ-2017-009/index.html
scrub_time_series = function(ts_file_bulk, confound_file_bulk, method="censor", n=NULL, n.nodes=NULL) {

  if (method=="censor") {
    badrows = lapply(confound_file_bulk, bad_ts)
    ts_file_bulk = lapply(ts_file_bulk, function(i) as.matrix(i))
    ts_file_bulk  = lapply(ts_file_bulk , function(l) l[,seq(1:n.nodes)])
    badrows.pre = lapply(badrows, function(l) l-1)
    badrows.pre = lapply(1:n, function(i) {badrows.pre[[i]][which(badrows.pre[[i]] < 1)] <- NA; badrows.pre[[i]]})
    badrows.pre = lapply(badrows.pre, na.omit)
    badrows.post = lapply(badrows, function(l) c(l+1, l+2))
    badrows.post = lapply(1:n, function(i) {badrows.post[[i]][which(badrows.post[[i]] > nrow(ts_file_bulk[[i]]))] <- NA; badrows.post[[i]]})
    badrows.post = lapply(badrows.post, na.omit)
    badrows = lapply(1:n, function(s) c(badrows[[s]], badrows.pre[[s]], badrows.post[[s]]))
    print("Removing Bad Time Points along with -1 through +2 time points.")
    badrows = lapply(badrows, function(i) as.integer(i))
    scrubbed_ts_bulk = lapply(ts_file_bulk, function(i) as.matrix(i))
    scrubbed_ts_bulk = pbapply::pblapply(1:n, function(i) {scrubbed_ts_bulk[[i]][badrows[[i]],1:n.nodes] <- NA; scrubbed_ts_bulk[[i]]})
    scrubbed_ts_bulk = lapply(scrubbed_ts_bulk, na.omit)
    scrubbed_ts_bulk = lapply(scrubbed_ts_bulk, function(i) as.matrix(i))
    return(scrubbed_ts_bulk)
    } else if (method=="interpolate") {
    badrows = lapply(confound_file_bulk, bad_ts)
    scrubbed_ts_bulk = lapply(ts_file_bulk, function(i) as.matrix(i))
    scrubbed_ts_bulk = lapply(scrubbed_ts_bulk, function(l) l[,seq(1:n.nodes)])
    badrows.pre = lapply(badrows, function(l) l-1)
    badrows.pre = lapply(1:n, function(i) {badrows.pre[[i]][which(badrows.pre[[i]] < 1)] <- NA; badrows.pre[[i]]})
    badrows.pre = lapply(badrows.pre, na.omit)
    badrows.post = lapply(badrows, function(l) c(l+1, l+2))
    badrows.post = lapply(1:n, function(i) {badrows.post[[i]][which(badrows.post[[i]] > nrow(ts_file_bulk[[i]]))] <- NA; badrows.post[[i]]})
    badrows.post = lapply(badrows.post, na.omit)
    badrows= lapply(1:n, function(s) c(badrows[[s]], badrows.pre[[s]], badrows.post[[s]]))
    rm(badrows.pre, badrows.post)

    print("Setting Bad Time Points along with -1 through +2 time points to NA.")
    scrubbed_ts_bulk = pbapply::pblapply(1:n, function(i) {scrubbed_ts_bulk[[i]][badrows[[i]],1:n.nodes] <- NA; scrubbed_ts_bulk[[i]]})

    print("Interpolating Missing Values")
    scrubbed_ts_bulk = pbapply::pblapply(scrubbed_ts_bulk, function(i) apply(i, 2, function(j) imputeTS::na.interpolation(j,option = "spline")))
    scrubbed_ts_bulk = lapply(scrubbed_ts_bulk, function(i) apply(i, 1, function(j) jitter(j, amount=.0001))) # adds a small amount of jitter row-wise to avoid interpolated time points being too identical across nodes
    scrubbed_ts_bulk = lapply(scrubbed_ts_bulk, function(i) apply(i, 2, function(j) jitter(j, amount=.0001))) # adds a small amount of jitter to avoid any interpolated columns being perfectly identical
    scrubbed_ts_bulk = lapply(scrubbed_ts_bulk, function(i) as.matrix(t(i)))
    return(scrubbed_ts_bulk)
    }

}
