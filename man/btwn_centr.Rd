% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/btwn_centr.R
\name{btwn_centr}
\alias{btwn_centr}
\title{Calculate betweenness centrality for a single graph.}
\usage{
btwn_centr(graph, weights = abs(E(graph)$weight))
}
\arguments{
\item{graph}{A network as an igraph object or matrix.}

\item{weights}{a vector of edge weights. Defaults to the absolute value of the input network.}
}
\value{
A vector of centrality scores for each node.
}
\description{
This function is a convenience wrapper for igraph's betweenness centrality function for consistency with rsfcNet's naming conventions.
}
\details{
Betweenness centrality gives the number of times a given node lies in the shortest path between two other nodes.
Betweenness centrality may be less natural to interpret in rsfc networks (see Power et al 2013)
but remains a popular centrality metric. See for example Wang, Zuo, & He (2010).
It is calculated with the formula below:

\eqn{b_i = \sum_{s \neq v \neq t}\frac{\sigma_{st}(i)}{\sigma_{st}}}
}
\examples{
betweenness = btwn_centr(graph)

}
\references{
Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4

Power, J. D., Schlaggar, B. L., Lessov-Schlaggar, C. N., & Petersen, S. E. (2013). Evidence for hubs in human functional brain networks. Neuron, 79(4), 10.1016/j.neuron.2013.07.035. http://doi.org/10.1016/j.neuron.2013.07.035

Wang, J., Zuo, X., & He, Y. (2010). Graph-Based Network Analysis of Resting-State Functional MRI. Frontiers in Systems Neuroscience, 4, 16. http://doi.org/10.3389/fnsys.2010.00016
}
\seealso{
\code{\link[rsfcNet]{btwn_centr_mult}}
}
\author{
Brandon Vaughan
}
