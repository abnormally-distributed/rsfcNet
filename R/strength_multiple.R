#'  Calculate signed strength statistics on a multiple graphs.
#'
#' This calculates strength statistics on a list of graphs.
#' @param graphs a list of networks in matrix format or as an igraph objects.
#' @param which defaults to the weighted strength metric "star". Can also be "pos", "neg", or "all".
#' @param col.names The names of each column (node labels)
#' @param row.names The names of each row (subject)
#' @param Scale Should Strength* be scaled to the range -1 to 1 (assuming correlations are the edge weights)? Defaults to "TRUE"
#' @return A matrix of the statistic for the chosen strength measure, or a list of matrices if which="all"
#' @export
#' @author Brandon Vaughan
#'
#' @details Strength star is the weighted combination of the positive connections and negative connections
#' in a network proposed by Rubinov and Sporns (2011). The statistic is calculated by the following formula:
#'
#'
#' \eqn{s_i^{*}=s_i^+-\left(\frac{s_i^-}{s_i^++s_i^-}\right)s_i^-}
#'
#'
#' where the positive and negative strength are respectively the sum of positive/negative weights:
#'
#'
#' \eqn{s_i^\pm = \sum_{j\neq i}^N w_{ij}}
#'
#'
#' When normalized to the -1 to 1 interval with scale=TRUE the positive and negative strength
#' are normalized with the following formula first then plugged into the above formula:
#'
#'
#' \eqn{s_i^\pm = \frac1{n-1}s_i^\pm}
#'
#'
#' @examples
#' strength = strength_multiple(graphs, which="star")
#'
#' @seealso
#' \code{\link[rsfcNet]{strength_signed}}
#' \code{\link[rsfcNet]{degree_centr}}
#' \code{\link[igraph]{strength}}
#' \code{\link[igraph]{degree}}
#'
#'
#' @references
#'
#' Fornito, A., Zalesky, A., & Bullmore, E. (2016). Node Degree and Strength. Chapter 4. Fundamentals of Brain Network Analysis, 115-136. doi:10.1016/B978-0-12-407908-3.00004-2
#'
#' Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. NeuroImage, 56(4), 2068-2079. doi:10.1016/j.neuroimage.2011.03.069
#`
strength_multiple = function(graphs, which="star", Scale=TRUE, col.names=NULL, row.names=NULL) {
  if (all(sapply(graphs, igraph::is.igraph)==TRUE)==TRUE) {

  graphs = lapply(graphs, function (g) as.matrix(igraph::as_adjacency_matrix(g, edges = FALSE, attr = "weight", sparse=TRUE)))

  }

  if (which=="all"){
    strength.stats = pbapply::pblapply(graphs, strength_signed)
    strength_pos =  sapply(strength.stats, function(i) i$positive_strength)
    strength_neg =  sapply(strength.stats, function(i) i$negative_strength)
    strength_signed = sapply(strength.stats, function(i) i$strength_star)
    strength_pos = data.frame(t(strength_pos))
    strength_neg = data.frame(t(strength_neg))
    strength_signed = data.frame(t(strength_signed))
    colnames(strength_pos) = col.names
    rownames(strength_pos) = row.names
    colnames(strength_neg) = col.names
    rownames(strength_neg) = row.names
    colnames(strength_signed) = col.names
    rownames(strength_signed) = row.names
    strength_statistics = list(strength_pos, strength_neg, strength_signed)
    return(strength_statistics)
  } else if (which=="star") {
    strength.star = t(pbapply::pbsapply(graphs, function(i) strength_signed(i, scale=Scale)$strength_star))
    colnames(strength.star) = col.names
    rownames(strength.star) = row.names
    return(strength.star)
  } else if (which=="pos") {
    strength.pos = t(pbapply::pbsapply(graphs, function(i) strength_signed(i, scale=Scale)$positive_strength))
    colnames(strength.pos) = col.names
    rownames(strength.pos) = row.names
    return(strength.pos)
  } else if (which=="neg") {
    strength.neg= t(pbapply::pbsapply(graphs, function(i) strength_signed(i, scale=Scale)$negative_strength))
    colnames(strength.neg) = col.names
    rownames(strength.neg) = row.names
    return(strength.neg)
  }

}
