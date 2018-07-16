#' Calculate Laplacian centrality for a single graph
#'
#' This calculates the Laplacian centrality metric for a graph.
#' @param graph A network as an igraph object. This function works much better on a network with thresholding applied.
#' @param prog.bar Should a progress bar be displayed? Defaults to TRUE.
#' @return A vector of the Laplacian centrality scores.
#' @author Brandon Vaughan
#' @references
#'
#' Pauls, S.D., & Remondini, D. (2012). A measure of centrality based on the spectrum of the Laplacian. Physical review. E, Statistical, nonlinear, and soft matter physics, 85 6 Pt 2, 066127.
#'
#' Qi, Xingqin, et al. (2012). Laplacian centrality: A new centrality measure for weighted networks. Information Sciences 194: 240-253.
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
#' laplace_centr(graph)
#'
#'  @seealso
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#' \code{\link[rsfcNet]{leverage_centr}}
#' \code{\link[rsfcNet]{leverage_centr_mult}}
#'
#' @export

laplace_centr <- function (graph, prog.bar=TRUE){

  ## This is just a function within a function to calculate unweighted strength star without the fluff of the
  ## function available for use in the package. It will return positive strength for all positive networks and
  ## returns degree for binary networks.

  strength_signed_special = function(graph){

    graph = as.matrix(igraph::as_adjacency_matrix(graph, edges = FALSE, attr = "weight", sparse=TRUE))
    m2 = graph
    diag(m2) = 0
    n.nodes=ncol(m2)
    pos.matrix <- m2
    pos.matrix[which(m2 < 0)] <- 0
    neg.matrix <- m2
    neg.matrix[which(m2 > 0)] <- 0
    pos.strength <- rowSums(pos.matrix)
    neg.strength <- abs(rowSums(neg.matrix))
    norm.strength <- cbind(pos.strength, neg.strength)
    strength.star <- as.vector(sapply(1:n.nodes, function(r) norm.strength[r,1] - (norm.strength[r,2]/(norm.strength[r,1]+norm.strength[r,2])*norm.strength[r,2])))
    strength.star[which(strength.star < 0)] <- 0
    strength.star

  }


  if (all(sign(E(graph)$weight)>=0)==TRUE) {      #If all edges are of weight greater than zero, igraph will be used to calculate strength/degree.
    nodes = V(graph)
    pb   <- txtProgressBar(1, as.numeric(length(V(graph))), style=3)
    centr <- integer();
    for (Node_i in V(graph)[nodes]) {
      if (prog.bar==TRUE) {
        setTxtProgressBar(pb, Node_i) }
      Neighbors <- neighborhood(graph, 1, nodes=Node_i)[[1]][-1]
      str <- 0
      for (nb in Neighbors) {
        str <- sum(str,igraph::strength(graph, nb, mode="all", loops = FALSE))
      }
      iStr <- igraph::strength(graph, Node_i, mode="all", loops = FALSE)
      centr <- append(centr, iStr^2 + iStr + 2*str )
    }
    return(centr)

} else if (all(sign(E(graph)$weight)>=0)==FALSE) {  #If any edges are negative, strength star will be used to calculate strength. Any negative strengths will be set to zero.

  nodes = V(graph)
  pb   <- txtProgressBar(1, as.numeric(length(V(graph))), style=3)
  centr <- integer();
  for (Node_i in V(graph)[nodes]) {
    if (prog.bar==TRUE) {
      setTxtProgressBar(pb, Node_i) }
    Neighbors <- neighborhood(graph, 1, nodes=Node_i)[[1]][-1]
    str <- 0
    for (nb in Neighbors) {
      str <- sum(str,strength_signed_special(graph)[nb])
    }
    iStr <- strength_signed_special(graph)[Node_i]
    centr <- append(centr, iStr^2 + iStr + 2*str )
  }
   return(centr)
 }
}