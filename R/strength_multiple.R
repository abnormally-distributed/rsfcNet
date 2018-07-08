#' Calculate node connection strength for weighted undirected networks (multiple subjects).
#'
#' This calculates the strength star metric proposed by Rubinov and Sporns (2011). Strength star is a weighted combination of the positive connections and negative connections in a correlation matrix. This function can also return positive and negative strength.
#' @param cormats a list of connectivity matrices.
#' @param which which = "all" returns the strength star metric from Rubinov and Sporns (2011) along with the positive strengths and negative strengths for each node. which="star", "pos", and "neg" return only the weighted strength, positive strength, and negative strength respectively.
#' @param col.names The names of each column (node labels)
#' @param row.names The names of each row (subject)
#' @param Scale Should Strength* be scaled to the range -1 to 1 (assuming correlations are the edge weights)? Defaults to Scale="TRUE"
#' @return A matrix of the statistic for the chosen strength measure, or a list of matrices if which="all"
#' @export
#' @author Brandon Vaughan
#' @examples
#' strength = strength_multiple(cormats, which="star")
#' @references
#' Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. NeuroImage, 56(4), 2068-2079. doi:10.1016/j.neuroimage.2011.03.069

strength_multiple = function(cormats, which="all", Scale="TRUE", col.names=.GlobalEnv$colnames, row.names=.GlobalEnv$rownames) {
  if (which=="all"){
    strength.stats = pbapply::pblapply(cormats, strength_star)
    strength_pos =  sapply(strength.stats, function(i) i$positive_strength)
    strength_neg =  sapply(strength.stats, function(i) i$negative_strength)
    strength_star = sapply(strength.stats, function(i) i$strength_star)
    strength_pos = data.frame(t(strength_pos))
    strength_neg = data.frame(t(strength_neg))
    strength_star = data.frame(t(strength_star))
    colnames(strength_pos) = col.names
    rownames(strength_pos) = row.names
    colnames(strength_neg) = col.names
    rownames(strength_neg) = row.names
    colnames(strength_star) = col.names
    rownames(strength_star) = row.names
    strength_statistics = list(strength_pos, strength_neg, strength_star)
    return(strength_statistics)
  } else if (which=="star") {
    strength.star = t(pbapply::pbsapply(cormats, function(i) strength_star(i, scale=Scale)$strength_star))
    colnames(strength.star) = col.names
    rownames(strength.star) = row.names
    return(strength.star)
  } else if (which=="pos") {
    strength.pos = t(pbapply::pbsapply(cormats, function(i) strength_star(i, scale=Scale)$positive_strength))
    colnames(strength.pos) = col.names
    rownames(strength.pos) = row.names
    return(strength.pos)
  } else if (which=="neg") {
    strength.neg= t(pbapply::pbsapply(cormats, function(i) strength_star(i, scale=Scale)$negative_strength))
    colnames(strength.neg) = col.names
    rownames(strength.neg) = row.names
    return(strength.neg)
  }

}
