#' Find current-flow closeness centrality for a single graph
#'
#' A variant of closeness centrality based on circuits
#' @details
#' Current flow centrality was developed based on the properties of electrical circuits. It is a variant of closeness
#' centrality that imagines edges resistors and the nodes as relays between resistors. Weighted edge information
#' is taken into account as the absolute value. Non-zero edge are conceptually "conductors" while zero edges are
#' "resistors."
#'
#' Current flow centrality is somewhat mathematically involved, but the final formula is
#'
#' \eqn{\displaystyle{\sum_j A_{i,j} (v_{i}^{(s,t)} - v_{j}^{(s,t)}) = u_{i}^{(s,t)}}$}.
#'
#' See the website in the references for a full explanation. The method involves calculating the laplacian representation
#' of a graph. This code is modified from the code in the centiserve package.
#'
#' @param graph An igraph object or correlation matrix
#' @return A numeric vector contaning the centrality scores for the selected vertices.
#' @author Brandon Vaughan
#' @author Mahdi Jalili
#'
#' @seealso
#' \code{\link[rsfcNet]{closeness_centr}}
#' \code{\link[rsfcNet]{laplace_centr}}
#'
#' @export
#' @references
#' Brandes U., Fleischer D. (2005) Centrality Measures Based on Current Flow. In: Diekert V., Durand B. (eds) STACS 2005. STACS 2005. Lecture Notes in Computer Science, vol 3404. Springer, Berlin, Heidelberg
#'
#' https://www.sci.unich.it/~francesc/teaching/network/flowcentrality.html
#'
#' @examples
#' current_centr(graph)
#'
current_centr = function (graph) {
  if (igraph::is.igraph(graph)==FALSE){
    graph = igraph::graph.adjacency(graph, mode="undirected", weighted=T, diag=FALSE)
  }
  vids = V(graph)
  E(graph)$weight = abs(E(graph)$weight)
  numVertices <- vcount(graph)
  resultVector <- double(numVertices)
  solveTemporaryVector <- matrix(nrow=numVertices - 1, ncol=1)
  kroneckerVector <- matrix(nrow=numVertices - 1, ncol=1)
  L  = as.matrix(graph.laplacian(graph))
  L <- L[-1,-1];
  for (i in 1:(numVertices - 1)) {
    solveTemporaryVector[,1] <- 0
    kroneckerVector[,1] <- 0
    kroneckerVector[i,1] <- 1
    solveTemporaryVector <- solve(L, kroneckerVector)
    resultVector[1] <- resultVector[1] + solveTemporaryVector[i, 1]
    resultVector[i + 1] <-  resultVector[i + 1] + solveTemporaryVector[i, 1]
    for (j in 1:(numVertices - 1)) {
      resultVector[i + 1] <- resultVector[i + 1] + (solveTemporaryVector[i, 1] - 2 * solveTemporaryVector[j, 1])
      resultVector[j + 1] <- resultVector[j + 1] + solveTemporaryVector[i, 1]
    }
  }
  resultVector[vids]
}


#' Find current-flow closeness centrality for a list of graphs
#'
#' A variant of closeness centrality based on circuits
#' @details
#' Current flow centrality was developed based on the properties of electrical circuits. It is a variant of closeness
#' centrality that imagines edges resistors and the nodes as relays between resistors. Weighted edge information
#' is taken into account as the absolute value. Non-zero edge are conceptually "conductors" while zero edges are
#' "resistors."
#'
#' Current flow centrality is somewhat mathematically involved, but the final formula is
#'
#' \eqn{\displaystyle{\sum_j A_{i,j} (v_{i}^{(s,t)} - v_{j}^{(s,t)}) = u_{i}^{(s,t)}}$}.
#'
#' See the website in the references for a full explanation. The method involves calculating the laplacian representation
#' of a graph. This code is modified from the code in the centiserve package.
#'
#' @param graphs A list of igraph objects or connectivity matrices
#' @param col.names The names of each column (node labels).
#' @param row.names The names of each row (subject).
#' @param parallel Should multiple cores be used? Defaults to FALSE. If TRUE, progress bar is not displayed. This is normal.
#' @param cores How many cores should be used? Defaults to recommended 1 less than number of CPU cores.

#'
#' @return A matrix contaning the centrality scores.
#' @author Brandon Vaughan
#'
#' @seealso
#' \code{\link[rsfcNet]{closeness_centr}}
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#' \code{\link[rsfcNet]{current_centr}}
#'
#' @export
#' @references
#' Brandes U., Fleischer D. (2005) Centrality Measures Based on Current Flow. In: Diekert V., Durand B. (eds) STACS 2005. STACS 2005. Lecture Notes in Computer Science, vol 3404. Springer, Berlin, Heidelberg
#'
#' https://www.sci.unich.it/~francesc/teaching/network/flowcentrality.html
#'
#' @examples
#' current_centr_mult(graph)
#'
current_centr_mult = function(graphs, col.names = NULL, row.names = NULL, parallel=TRUE, cores=NA){
  if (parallel==FALSE){
  current.centrality = pbapply::pbsapply(graphs, function(x) current_centr(x))
  current.centrality = t(current.centrality)
  colnames(current.centrality) = col.names
  rownames(current.centrality) = row.names
  return(current.centrality)
  } else if (parallel==TRUE) {
    if (is.na(cores)==TRUE) {
      cl = as.integer(parallel::detectCores()-1)
      cl = parallel::makeCluster(cl)
    } else if (is.na(cores)==FALSE) {
      cl = as.integer(cores)
      cl = parallel::makeCluster(cl)
    }

    clusterEvalQ(cl, {
      library(parallel)
      library(igraph)
    })

    current.centrality  = parallel::parSapply(cl, graphs, function(x) rsfcNet::current_centr(x))
    stopCluster(cl)
    current.centrality = t(current.centrality)
    colnames(current.centrality) = col.names
    rownames(current.centrality) = row.names
    return(current.centrality)
  }
}
