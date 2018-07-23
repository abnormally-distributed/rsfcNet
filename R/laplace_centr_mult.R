#' Calculate Laplacian centrality for a single graph
#'
#' This calculates the Laplacian centrality metric for a graph.
#' @param graphs a list of igraph objects.
#' @param parallel Should multiple cores be used? Defaults to FALSE. If TRUE, progress bar is not displayed. This is normal.
#' @param cores How many cores should be used? Defaults to recommended 1 less than number of CPU cores.
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @export
#' @author Brandon Vaughan
#'
#' @details
#' To first understand Laplacian centrality the concept of the Laplacian matrix must be understood.
#' The Laplacian matrix (or graph Laplacian) is the adjacency
#' (or correlation or other weighted connectivity) matrix subtracted from the
#' degree (or strength) matrix. In the strength matrix all entries are zero except for the diagonal,
#' which has as an entry the number corresponding to the degree or strength of node_i.
#'
#' Mathematically this is represented as
#'
#' \eqn{L(G) = X(G) - W(G)}
#'
#' or
#'
#' \eqn{L(G) = D(G) - A(G)}
#'
#' The energy of L(G) is used in physics to measure things such as the diffusion
#' of energy through a system and is shown in the formula below:
#'
#' \eqn{E_L(G)=\sum_{i=1}^nx_i^2+2\sum_{i<j}w_{i\text{,}j}^2}
#'
#' Laplacian centrality extends the idea by studying what happens to the ability of energy (or information)
#' to difuse through a network if a node is removed. The centrality score for each node is calculated by
#' the expected drop in energy for the graph Laplacian when node_i is removed. A method of calculating this
#' is given by the formula below:
#'
#' \eqn{\Delta E_{n}=s_{G}^{2}(n)+s_{G}(n)+2\sum_{n _{i}\in N(n)}s_{G}(n_{i})}
#'
#' This measure is very appropriate to the study of brain networks as it can charaterize the expected change
#' in network energy if a region where removed from the brain.
#'
#' @examples
#' laplace_centr_mult(graphs, parallel=FALSE) # If you only have two cores parallel=FALSE is a good idea.
#'
#' @seealso
#' \code{\link[rsfcNet]{laplace_centr}}
#' \code{\link[rsfcNet]{leverage_centr}}
#' \code{\link[rsfcNet]{leverage_centr_mult}}
#'
#' @references
#' Pauls, S.D., & Remondini, D. (2012). A measure of centrality based on the spectrum of the Laplacian. Physical review. E, Statistical, nonlinear, and soft matter physics, 85 6 Pt 2, 066127.
#'
#' Qi, Xingqin, et al. (2012). Laplacian centrality: A new centrality measure for weighted networks. Information Sciences 194: 240-253.
#'
laplace_centr_mult = function (graphs, col.names = NULL, row.names = NULL, parallel=FALSE, cores=NA)
    {
if (parallel==FALSE) {
  laplace_centr = pbapply::pbsapply(graphs, function(x) laplace_centr(x,prog.bar=FALSE))
  laplace_centr = t(laplace_centr)
  colnames(laplace_centr) = col.names
  rownames(laplace_centr) = row.names
  laplace_centr = as.matrix(laplace_centr)
  return(laplace_centr)
}  else if (parallel==TRUE) {

  if (is.na(cores)==TRUE) {
    cl = as.integer(parallel::detectCores()-1)
    cl = makeCluster(cl)
  } else if (is.na(cores)==FALSE) {
    cl = as.integer(cores)
    cl = makeCluster(cl)
  }

  clusterEvalQ(cl, {
    library(igraph)
  })

  laplace_centr = parallel::parSapply(cl, graphs, function(x) rsfcNet::laplace_centr(x,prog.bar=FALSE))
  stopCluster(cl)
  laplace_centr = t(laplace_centr)
  colnames(laplace_centr) = col.names
  rownames(laplace_centr) = row.names
  laplace_centr = as.matrix(laplace_centr)
  return(laplace_centr)
  }
}


