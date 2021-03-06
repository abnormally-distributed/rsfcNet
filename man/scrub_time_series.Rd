% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrub_time_series.R
\name{scrub_time_series}
\alias{scrub_time_series}
\title{"Scrubs" Bad Time Points from Time Series Matrix}
\usage{
scrub_time_series(
  ts_file_bulk,
  confound_file_bulk,
  method = "censor",
  n = NULL,
  n.nodes = NULL
)
}
\arguments{
\item{ts_file_bulk}{list of time series matrices.}

\item{confound_file_bulk}{list of confound files for each subject}

\item{method}{method="censor" uses censoring/scrubbing while interpolation fills in missing time points column-wise with smooth polynomials.}

\item{n}{The The number of subjects in the list.}

\item{n.nodes}{The The number of nodes in the network.}
}
\value{
The time series matrix for each subject.
}
\description{
This function takes a list of time series files and a list of confound matrices and returns the scrubbed time series files. The default method is "censor". Censoring takes the preceding time point and the following two time points of each bad time point in the confound matrices and deletes them as described in Power et al, 2012. If there are relatively few bad time points in the subjects method="interpolate" may be tried.
}
\examples{
scrubbed_ts = scrub_time_series(ts, confounds)
}
\references{
Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage, 59(3), 2142–2154. http://doi.org/10.1016/j.neuroimage.2011.10.018

Moritz S, Bartz-Beielstein T (2017). “imputeTS: Time Series Missing Value Imputation in R.” The R Journal, 9(1), 207–218. https://journal.r-project.org/archive/2017/RJ-2017-009/index.html
}
\author{
Brandon Vaughan
}
