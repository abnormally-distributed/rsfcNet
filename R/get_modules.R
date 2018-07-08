#' Get the modularity and module partitions of a weighted undirected network.
#'
#' This function takes an igraph object and uses igraph to calculate the modularity and community structure of the graph.
#' @param graph an igraph object
#' @param method Can be "fast_greedy", "louvain", "walktrap", or "eigenvector".
#' @param step.length The walk length for the walktrap algorithm. Best to set it to the number of modules expected on a theoretical basis.
#' @return A communities (igraph) object containing a membership vector and modularity statistic.
#' @export
#' @author Brandon Vaughan
#' @details
#' This function simplifies the process over using igraph directly. All edges are converted to their absolute value. Only certain algorithms are allowed for computational effeciency when dealing with large graphs. For method="fast_greedy" each node initially is its own module, and the algorithm combines modules iteratively. At each step the modules are merged with the one which would give the largest increase in the modularity score. The merging stops when modularity can be increased no further. It is considered a "bottom-up" approach (Clauset, Newman, & Moore, 2004). The "fast_greedy" method has a tendency to under-estimate the number of modules when the module size is small (Yang, Algesheimer, & Tessone, 2016). method="louvain" works in a similar manner by beginning with each node assigned to its own module and merging them. After each merge each module is treated as a single node until it reaches a point where modularity is maximized (Blondel et al, 2008). By contrast, method="eigenvector" uses a "top-down" approach that begins with the entire network as one module and iteratively tests different partitions until the modularity score can be increased no further (Newman, 2006). This process can be thought of as analagous to a principal components analysis where the component on which a given node has the strongest loading defines the module membership. igraph also offers Spinglass and Edge Betweenness algorithms, but the computational time increases dramatically with network size and thus are not recommended for larger networks (Yang, Algesheimer, & Tessone, 2016). Although a typical atlas-based network doesn't qualify as a "large network" in the typical sense, run time for Spinglass and Edge-Betweenness increases dramatically fast such that even a 264 node atlas is too large. Gates et al (2016) conclude that for applications to neuroscience the louvain and walktrap methods are both effecient and reasonably accurate. The fast greedy method performs well in the presence of a sparse connectivity matrix (which shrinkage estimated partial correlation matrices often are) so long as the number of nodes is more than 100 (Gates et al, 2016). For more information see Chapter 9 of Fornito et al (2016).
#' @examples
#' modules = get_modules(graph, method="louvain")
#'
#' #Applied to a list of graphs:
#' module_list = lapply(graphs, function(i) get_modules(i, method="louvain"))
#'
#' @references
#' Blondel, V. D., Guillaume, J., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. Journal of Statistical Mechanics: Theory and Experiment, 2008(10). doi:10.1088/1742-5468/2008/10/p10008
#'
#' Clauset, A., Newman, M. E., & Moore, C. (2004). Finding community structure in very large networks. Physical Review E, 70(6). doi:10.1103/physreve.70.066111
#'
#` Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 9. Fundamentals of Brain Network Analysis, 303-354. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Gates, K. M., Henry, T., Steinley, D., & Fair, D. A. (2016). A Monte Carlo Evaluation of Weighted Community Detection Algorithms. Frontiers in Neuroinformatics, 10. doi:10.3389/fninf.2016.00045
#'
#' Newman, M. E. (2006). Finding community structure in networks using the eigenvectors of matrices. Physical Review E, 74(3). doi:10.1103/physreve.74.036104
#'
#' Yang, Z., Algesheimer, R., & Tessone, C. J. (2016). A Comparative Analysis of Community Detection Algorithms on Artificial Networks. Scientific Reports, 6(1). doi:10.1038/srep30750
#'
get_modules = function(graph, method="fast_greedy", step.length) {
    if (method=="fast_greedy") {
    greedgraph = cluster_fast_greedy(graph, weights=abs(E(graph)$weight), modularity = TRUE, membership = TRUE)
    greedgraph$modularity = modularity(greedgraph) #seems redundant with modularity=TRUE but otherwise it returns a vector of the modularity after each merge in the algorithm.
    return(greedgraph)
  } else if (method=="eigenvector") {
    eigengraph = cluster_leading_eigen(graph, weight=abs(E(graph)$weight))
    eigengraph$modularity = modularity(eigengraph)
    eigengraph$membership = membership(eigengraph)
    return(eigengraph)
  } else if (method=="louvain") {
    louvgraph = cluster_louvain(graph, weights =abs(E(graph)$weight))
    louvgraph$modularity = modularity(louvgraph)
    louvgraph$membership = membership(louvgraph)
    return(louvgraph)
  } else if (method=="walktrap") {
    walkgraph = cluster_walktrap(graph, weights =abs(E(graph)$weight), steps=step.length)
    walkgraph$modularity = modularity(walkgraph)
    walkgraph$membership = membership(walkgraph)
    return(walkgraph)
  }
}
