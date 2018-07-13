#' Calculate signed strength statistics on a single graph.
#'
#' This calculates strength statistics on a single graph.
#' @param matrix A network in matrix format or as an igraph object.
#' @param scale Defaults to TRUE to return the scaled statistic. Can be set to FALSE.
#' @return The weighted strength metric along with the positive and negative strength.
#' @export
#' @author Brandon Vaughan
#'
#' @details Strength star is the weighted combination of the positive connections and negative connections
#' in a network proposed by Rubinov and Sporns (2011). The statistic is calculated by the following formula:
#' \eqn{s_i^{*}=s_i^+-\left(\frac{s_i^-}{s_i^++s_i^-}\right)s_i^-}
#' When normalized to the -1 to 1 interval with scale=TRUE the positive and negative strength
#' are normalized with the following formula first then plugged into the above formula:
#' \eqn{s_i^\pm = \frac1{n-1}s_i^\pm}
#'
#'
#' @examples strength_star(sub001_cormat, scale=TRUE)
#' @references Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. NeuroImage, 56(4), 2068-2079. doi:10.1016/j.neuroimage.2011.03.069
strength_signed = function(matrix, scale=TRUE){
  if (is.igraph(matrix)=="TRUE"){
  matrix = as.matrix(igraph::as_adjacency_matrix(matrix, edges = FALSE, attr = "weight", sparse=TRUE)) }
  m2 = matrix
  diag(m2) = 0
  n.nodes=ncol(m2)
  pos.matrix <- m2
  pos.matrix[which(m2 < 0)] <- 0
  neg.matrix <- m2
  neg.matrix[which(m2 > 0)] <- 0
  pos.strength <- rowSums(pos.matrix)
  neg.strength <- abs(rowSums(neg.matrix))

  if (scale==TRUE) {
    norm.pos<- (1/(n.nodes-1))*pos.strength
    norm.neg<- (1/(n.nodes-1))*neg.strength
    norm.strength <- cbind(norm.pos, norm.neg)
    strength.star <- as.vector(sapply(1:n.nodes, function(r) norm.strength[r,1] - (norm.strength[r,2]/(norm.strength[r,1]+norm.strength[r,2])*norm.strength[r,2])))
    return(cbind.data.frame(positive_strength=norm.pos , negative_strength=norm.neg, strength_star = strength.star))
  } else if (scale==FALSE) {
    norm.strength <- cbind(pos.strength, neg.strength)
    strength.star <- as.vector(sapply(1:n.nodes, function(r) norm.strength[r,1] - (norm.strength[r,2]/(norm.strength[r,1]+norm.strength[r,2])*norm.strength[r,2])))
    return(cbind.data.frame(positive_strength=pos.strength , negative_strength=neg.strength, strength_star = strength.star))
  }
}


