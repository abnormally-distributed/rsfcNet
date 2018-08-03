#' Calculate several module/community based network statistics
#'
#' This calculates the (within-module) strength star and diversity star metrics proposed by Rubinov and Sporns (2011) along with several other metrics.
#' @param graph a single network in matrix format or as an igraph object
#' @param modules a communities object from an igraph clustering function
#' @param n.nodes The The number of nodes in each network.
#' @return A matrix containing several modularity/community structure based centrality measures.
#' @export
#' @author Brandon Vaughan
#'
#' @details This calculates the within module strength (weighted degree) z-score (using the weighted
#' signed "strength star" metric) and "diversity star" metric proposed by Rubinov and Sporns (2011).
#' This also returns the participation coefficient weighted by positive and negative connections in the
#' same manner as diversity. Also returned are the percentage of positive/negative weighted degree
#' a node has within-module out of its total weighted degree.

#' A node that is well-integrated into a community/module should arguably have
#' a high percentage of its total positive connections inside the module,
#' and a relatively low percentage of its negative connections within the module.
#'
#'
#' The formula used for diversity:
#'
#' \eqn{h_i = -\frac{1}{log(m)}  \sum_{u=1}^{N_M} \Bigg( \left ( \frac{s_{iu}}{s_i} \right ) \cdot log  \left ( \frac{s_{iu}}{s_i} \right ) \Bigg)}
#'
#' The formula used for participation:
#'
#' \eqn{p_i = 1 - \sum_{s=1}^{N_M} \left ( \frac{s_{iu}}{s_i} \right )^2}
#'
#'
#' The weighting for diversity is calculated as follows (and participation the same formula):
#'
#' \eqn{h_{i}^{*} = h_{i}^{+} - \Bigg( \frac{s_i^{-}}{s_i^{+}+s_i^{-}}\Bigg) h_{i}^{-}}
#'
#' @examples
#' \dontrun{
#' module_stats = module_connectivity(graph, module)
#' }
#'
#' @references
#'
#' Power, J. D., Schlaggar, B. L., Lessov-Schlaggar, C. N., & Petersen, S. E. (2013). Evidence for hubs in human functional brain networks. Neuron, 79(4), 10.1016/j.neuron.2013.07.035
#'
#' Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. NeuroImage, 56(4), 2068-2079. doi:10.1016/j.neuroimage.2011.03.069
#'

module_connectivity = function(graph, modules, n.nodes=NULL, consensus.output=FALSE) {

  if (is.igraph(graph)=="TRUE"){
    graph = as.matrix(igraph::as_adjacency_matrix(graph, edges = FALSE, attr = "weight", sparse=TRUE)) }

  matrix.mod = graph
  diag(matrix.mod) <- 0

  if (consensus.output==TRUE) {
    module_citizens = modules$groups
    module_matrices = lapply(1:length(unique(module_citizens)), function(m) matrix.mod[module_citizens[[m]],module_citizens[[m]]])
  } else {
    module_citizens = lapply(1:length(unique(modules$membership)), function(u) which(modules$membership==u))
    module_matrices = lapply(1:length(unique(modules$membership)), function(m) matrix.mod[module_citizens[[m]],module_citizens[[m]]])
  }

  strength_modules = lapply(module_matrices, function(m) strength_signed(m, scale=FALSE))
  positive_strengths = lapply(strength_modules, function(m) m$positive_strength)
  negative_strengths = lapply(strength_modules, function(m) m$negative_strength)
  negative_strengths = lapply(negative_strengths, function(n) {n[n == "NaN"] <-  0; n})
  module_strength_star =  lapply(strength_modules, function(m) m$strength_star)

  # Sort Nodes into original order
  sorted = cbind(module.pos= unlist(positive_strengths), module.neg= unlist(negative_strengths), module.star=unlist(module_strength_star), Node.ID = unlist(module_citizens))
  sorted = sorted[order(sorted[,4]),]

  #Get Node Strengths Across Network
  strength_star = strength_signed(graph, scale = FALSE)
  sorted.net = cbind.data.frame(sorted, pos.network = strength_star$positive_strength, neg.network = strength_star$negative_strength, star.network = strength_star$strength_star)

  hplus = sorted.net$module.pos/sorted.net$pos.network
  hplus = (-1 * 1/log(length(modules))) * (hplus*log(hplus))

  hneg= sorted.net$module.neg/sorted.net$neg.network
  hneg = (-1 * 1/log(length(modules))) * (hplus*log(hneg))
  hneg[hneg == Inf] <-  0
  hneg[hneg == -Inf] <- 0

  if (all(graph >= 0)=="TRUE") {
    diversity = hplus
  } else {
    diversity = hplus - ((sorted.net$neg.network/(sorted.net$neg.network+sorted.net$pos.network))*hneg)
  }

  ## Get Participation Coefficient

  part_pos = (1-(sorted.net$module.pos/sorted.net$pos.network))^2
  part_neg = (1-(sorted.net$module.neg/sorted.net$neg.network))^2

  if (all(graph >= 0)=="TRUE") {
    participation = part_pos
  } else {
    participation = part_pos - ((sorted.net$neg.network/(sorted.net$neg.network+sorted.net$pos.network))*part_neg)
  }

  PctTotalPosCon = sorted.net$module.pos/sorted.net$pos.network
  PctTotalNegCon = sorted.net$module.neg/sorted.net$neg.network
  participation.stats = cbind.data.frame(PctTotalPosCon, PctTotalNegCon, participation)

  ## Within Module Z-Score (based on strength star)

  zscores = lapply(module_strength_star, scale)
  zscores = cbind(unlist(zscores), Node.ID = unlist(module_citizens), unlist(module_strength_star))
  zscores = zscores[order(zscores[,2]),]

  module.statistics = cbind(within.module.z.score=zscores[,1],participation.stats, diversity=diversity)
  module.statistics = as.matrix(module.statistics)
  return(module.statistics)
}

