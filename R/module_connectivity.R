#' Calculate several module/community based network statistics
#'
#' This calculates the (within-module) strength star and diversity star metrics proposed by Rubinov and Sporns (2011) along with several other metrics.
#' @param graph a single network in matrix format or as an igraph object
#' @param modules a communities object from an igraph clustering function
#' @param n.nodes The The number of nodes in each network.
#' @param scale Whether or not to scale the strength of each node into a 0-1 range (TRUE) or simply return the weighted sum of connections (FALSE).
#' @return A matrix containing several modularity/community structure based centrality measures.
#' @export
#' @author Brandon Vaughan
#'
#' @details This calculates the (within-module) strength star and diversity star
#'  metrics proposed by Rubinov and Sporns (2011).
#' Also returns Within-Module Z-Score, Within-Module Strength, Diversity-Star,
#' Participation Coefficient, and percentage of total connections that exist within-module.
#' Z-Scores and Participation Coefficient measures are based on the strength (star) of the nodes.
#' In the event a binary network is fed, diversity will not be calculated.
#' Participation will be the more appropriate measure.
#' The percentage of positive intra-module connections for each
#' node out of the total positive-weight connetions
#' is also provided, along with the same measure for negative connections. A node that is well-integrated
#' into a community/module should arguably have more a high percentage of its total
#' positive connections inside the module,
#' and a relatively low percentage of its negative connections within the module.
#'
#' Diversity star is calculated by weighting the diversity coeffecient for positive and negative
#' connections.
#'
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

module_connectivity = function(graph, modules, scale=FALSE, n.nodes=NULL, consensus.output=FALSE) {

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

  strength_modules = lapply(module_matrices, function(m) strength_signed(m, scale=scale))
  positive_strengths = lapply(strength_modules, function(m) m$positive_strength)
  negative_strengths = lapply(strength_modules, function(m) m$negative_strength)
  negative_strengths = lapply(negative_strengths, function(n) {n[n == "NaN"] <-  0; n})
  module_strength_star =  lapply(strength_modules, function(m) m$strength_star)

  # Sort Nodes into original order
  sorted = cbind(module.pos= unlist(positive_strengths), module.neg= unlist(negative_strengths), module.star=unlist(module_strength_star), Node.ID = unlist(module_citizens))
  sorted = sorted[order(sorted[,4]),]

  #Get Node Strengths Across Network
  strength_star = strength_signed(graph, scale = scale)
  sorted.net = cbind.data.frame(sorted, pos.network = strength_star$positive_strength, neg.network = strength_star$negative_strength, star.network = strength_star$strength_star)

  hplus = sorted.net$module.pos/sorted.net$pos.network
  hplus = (-1 * 1/log(length(modules))) * (hplus*log(hplus))

  hneg= sorted.net$module.neg/sorted.net$neg.network
  hneg = (-1 * 1/log(length(modules))) * (hplus*log(hneg))
  hneg[hneg == Inf] <-  0
  hneg[hneg == -Inf] <- 0

  if (all(graph >= 0)=="TRUE") {
    diversity_star = diversity(graph)
  } else {
    diversity_star = hplus - (sorted.net$neg.network/(sorted.net$neg.network+sorted.net$pos.network))*hneg
  }

  ## Get Participation Coefficient

  participation = (1-(sorted.net$module.star/sorted.net$star.network))^2
  PctTotalPosCon = sorted.net$module.pos/sorted.net$pos.network
  PctTotalNegCon = sorted.net$module.neg/sorted.net$neg.network
  participation.stats = cbind.data.frame(module.strength=sorted.net$module.star, network.strength=sorted.net$star.network,PctTotalPosCon, PctTotalNegCon, participation)

  ## Within Module Z-Score (based on strength star)

  zscores = lapply(module_strength_star, scale)
  zscores = cbind(unlist(zscores), Node.ID = unlist(module_citizens), unlist(module_strength_star))
  zscores = zscores[order(zscores[,2]),]

  module.statistics = cbind(within.module.z.score=zscores[,1],participation.stats, diversity_star=diversity_star)
  module.statistics = as.matrix(module.statistics)
  return(module.statistics)
}

