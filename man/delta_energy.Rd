% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectral_measures.R
\name{delta_energy}
\alias{delta_energy}
\title{Calculate the nodewise-deletion change in graph energy for a single graph}
\usage{
delta_energy(graph)
}
\arguments{
\item{graph}{A network as an igraph object}
}
\value{
A vector of centrality scores for each node.
}
\description{
This function calculates the change in graph energy observed due to the deletion of each node in the graph.
}
\details{
Fornito, Zalesky, and Bullmore define a delta centrality measure as one that measures the change in
a global property of the graph that occurs due to the deletion of a node or edge (Fornito et al, 2016).
One informative global property is the energy of a graph. Graph energy was originally
applied in  organic chemistry to quantify the stability of molecular orbitals
associated with pi-electrons (Li, Shi, & Gutman, 2012). The graph energy informs about the connectivity of the graph
as a whole. A graph energy of zero means the nodes are not connected at all. Therefore,
to understand the importance of a particular node to the connectivity of the graph
a node can be removed and the change in graph energy used as a measure of node importance.
}
\examples{

delta_e = delta_energy(Graph)

}
\references{
Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4

Li, X.; Shi, Y.; Gutman, I. (2012), Graph Energy, New York: Springer, ISBN 978-1-4614-4219-6.
}
\seealso{
\code{\link[rsfcNet]{laplace_centr}}
\code{\link[rsfcNet]{graph_energy}}
\code{\link[rsfcNet]{fiedler_value}}
}
\author{
Brandon Vaughan
}
