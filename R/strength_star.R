#' Calculate node connection strength for weighted undirected networks.
#'
#' This calculates the strength star metric proposed by Rubinov and Sporns (2011). Strength star is a weighted combination of the positive connections and negative connections in a correlation matrix.
#' @param matrix a single connectivity matrix .
#' @param scale Whether or not to scale the strength of each node into a 0-1 range (TRUE) or simply return the weighted sum of connections (FALSE).
#' @return The scrubbed time series matrix.
#' @export
#' @author Brandon Vaughan
#' @examples strength_star(sub001_cormat, scale=TRUE)
#' @references Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. NeuroImage, 56(4), 2068-2079. doi:10.1016/j.neuroimage.2011.03.069
strength_star = function(matrix, scale=TRUE){
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


