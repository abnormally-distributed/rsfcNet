#' Calculate signed clustering coeffecient for a single graph
#'
#' Calculate the signed clustering coeffecient (local transitivity) for a single graph using
#' either the Zhang & Horvath method or the Constantini & Perugini method.
#'
#' @param graph A graph as an igraph object or matrix..
#' @param method Either "Constantini" (the default) or "Zhang"
#' @return The connectivity matrices.
#' @export
#' @author Giulio Costantini, Sacha Epskamp, Brandon Vaughan
#'
#' @details
#' Calculate the signed clustering coeffecient (local transitivity) for a single graph using
#' either the Zhang & Horvath (2005) method or the Constantini & Perugini (2014) method.
#' The local transitivity (or clustering coeffecient) of a node is the fraction of
#' edges a node forms with its neighbors out of the the number of edges it would
#' take to make complete triangles. This code was adapted from the code in the qgraph package
#' with modifications to make the output as in the matlab script
#' clustering_coef_w_sign in the Brain Connectivity Toolbox. The Onella method is excluded due to similarities with the Barrat (available as the local_clustcoef
#' function in this package) and Zhang methods.
#'
#' @seealso
#' \link[rsfcNet]{local_trans}
#' \link[igraph]{transitivity}
#' \link[rsfcNet]{clustcoef_signed_mult}
#'
#' @references
#' \href{https://sites.google.com/site/bctnet/measures/list}{Brain Connectivity Toolbox}
#'
#' Costantini, G., & Perugini, M. (2014). Generalization of Clustering Coefficients to Signed Correlation Networks. PLoS ONE, 9(2), e88669. http://doi.org/10.1371/journal.pone.0088669
#'
#' Zhang, B., & Horvath, S. (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology, 4(1).

clustcoef_signed = function(graph, method="Constantini")
{
  if (is.igraph(graph)=="TRUE") {
    W = as.matrix(as_adjacency_matrix(graph, attr="weight"))
    diag(W) <- 0
  } else if (is.igraph(graph)=="FALSE"){
    W = graph
    diag(W) <- 0
  }

  if (method == "Constantini") {
    a_W <- abs(W)
    num <- diag(W %*% W %*% W)
    a_num <- diag(a_W %*% a_W %*% a_W)
    den <- colSums(a_W)^2 - colSums(W^2)
    cZ <- num/den
    a_cZ <- a_num/den
    return(cZ)

  } else if (method=="Zhang") {
    W_pos = W
    W_pos[which(W_pos < 0)] = 0

    a_W_pos <- abs(W_pos)
    num_pos <- diag(W_pos %*% W_pos %*% W_pos)
    a_num_pos <- diag(a_W_pos %*% a_W_pos %*% a_W_pos)
    den_pos <- colSums(a_W_pos)^2 - colSums(W_pos^2)
    cZ_pos <- num_pos/den_pos
    a_cZ_pos <- a_num_pos/den_pos

    W_neg = W
    W_neg[which(W_neg > 0)] = 0
    W_neg = abs(W_neg)

    a_W_neg <- abs(W_neg)
    num_neg <- diag(W_neg %*% W_neg %*% W_neg)
    a_num_neg <- diag(a_W_neg %*% a_W_neg %*% a_W_neg)
    den_neg <- colSums(a_W_neg)^2 - colSums(W_neg^2)
    cZ_neg <- num_neg/den_neg
    a_cZ_neg <- a_num_neg/den_neg

    zhang = (a_cZ_neg + a_cZ_pos)/2
    return(zhang)
  }
}



