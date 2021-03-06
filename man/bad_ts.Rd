% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bad_ts.R
\name{bad_ts}
\alias{bad_ts}
\title{Read Bad Time Points from Confound Matrix}
\usage{
bad_ts(confound_matrix)
}
\arguments{
\item{confound_matrix}{The matrix of confounds.}
}
\value{
A vector indicating which time points are in need of scrubbing.
}
\description{
This identifies time points that must be scrubbed away. The function expects a confound matrix to be in the form provided by the fsl_motion_outliers function in the FSL software library.
}
\examples{
bad_ts(subj001_confoundm)
}
\references{
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLMotionOutliers
}
\author{
Brandon Vaughan
}
