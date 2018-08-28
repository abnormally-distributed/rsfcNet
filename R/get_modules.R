#' Get the modularity and module partitions of a weighted undirected network.
#'
#' This function takes an igraph object and uses igraph to calculate the modularity and community structure of the graph.
#' @param graph an igraph object
#' @param method Can be "fast_greedy", "louvain", "walktrap", "eigenvector", or "iterative_spinglass"
#' @param step.length The walk length for the walktrap algorithm.
#' @param upper.limit For method="iterative_spinglass", the upper limit for the number of modules expected.
#' @param signed For "iterative_spinglass", is the network signed?
#' @param repetitions Number of iterations for "iterative_spinglass"
#' @return A communities (igraph) object containing a membership vector and modularity statistic. For method="iterative_spinglass"
#' only a matrix of community label assignments for each run is returned.
#' @export
#' @author Brandon Vaughan
#' @details
#' All functions rely on the igraph library. All edges are converted to their absolute value for all but the signed iterative
#' spinglass method out of necessity. If you don't want this to happen, create versions of your graphs where the negative values are set to zero.
#'
#' In the "fast_greedy" method each node initially is its own module, and the algorithm combines modules iteratively.
#' At each step the modules are merged with the one which would give the largest increase in the modularity score.
#' The merging stops when modularity can be increased no further. It is considered a "bottom-up" approach
#' (Clauset, Newman, & Moore, 2004).
#' "fast_greedy" has a tendency to under-estimate the number of modules when modules are small in size
#' (Yang, Algesheimer, & Tessone, 2016).
#'
#'
#' The "louvain" methd works in a similar manner by beginning with each node assigned to its own module and merging them.
#' After each merge each module is treated as a single node until it reaches a point where modularity is maximized (Blondel et al, 2008).
#'
#'
#' By contrast, method="eigenvector" uses a "top-down" approach that begins with the entire network as one module and iteratively tests different partitions until the modularity score can be increased no further (Newman, 2006).
#' This process can be thought of as analagous to a principal components analysis where the component on which a given node has the strongest
#' loading defines the module membership.
#'
#'
#' The "walktrap" uses a very fast random walk algorithm algorithm to determine modular structure. Works well for dense weighted graphs and is
#' among the better algorithms (Gates et al., 2016; Yang et al., 2016). Highly recommend for applying to a list of graphs (see example).
#'
#'
#' The "iterative_spinglass" method uses a multi-iterative version of igraph's spinglass algorithm. A spinglass is
#' a material in which the north-south alignments of magnetic fields are disorganized. The community structure of
#' the network is interpreted as the spin configuration that minimizes the energy of the spin glass. Community indices
#' are given by the spin states. The ground state configuration gives a concise definition
#' of communities as cohesive subgroups (ie, groups with the same spin)
#' (Traag & Bruggeman, 2009; Wang et al., 2013; Zhang & Moore,2014) This has the advantage of working well with
#' signed networks. If your network is all-positive, however, setting signed=FALSE will allow a slightly
#' faster version of the algorithm to be run. If signed=FALSE and a signed network is given, the absolute
#' value of the network will be used. The option "repetitions" determines how many runs will be used. Bear in mind
#' that if the upper limit is more than 4 or so, there may be significant instability in the modules unless many
#' repetitions are used. Returned is a matrix of module assignments for uses with the modularity_consensus function.
#' This is a very slow method and will require a good deal of RAM and CPU power.
#'
#' @examples
#' modules = get_modules(graph, method="louvain")
#'
#' #Applied to a list of graphs:
#' module_list = lapply(graphs, function(i) get_modules(i, method="walktrap"))
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
#' Traag, V. & Bruggeman, J., (2009) Community detection in networks with positive and negative links. Physical Review E, 80(3) doi:10.1103/PhysRevE.80.036115.
#'
#' Wang, Z., Hu, Y., Xiao, W., & Ge, B. (2013). Overlapping community detection using a generative model for networks. Physica A: Statistical Mechanics and its Applications, 392(20), 5218-5230. doi:10.1016/j.physa.2013.06.038
#'
#' Yang, Z., Algesheimer, R., & Tessone, C. J. (2016). A Comparative Analysis of Community Detection Algorithms on Artificial Networks. Scientific Reports, 6(1). doi:10.1038/srep30750
#'
#' Zhang, P., Moore, C. (2014) Scalable detection of statistically significant communities and hierarchies, using message passing for modularity. Proceedings of the National Academy of Sciences Dec 2014, 111 (51) 18144-18149; DOI: 10.1073/pnas.1409770111
#'
get_modules = function(graph, method="fast_greedy", upper.limit=7, step.length=4, signed=TRUE, repetitions=25) {
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
    walkgraph = cluster_walktrap(graph, weights =abs(E(graph)$weight), steps=upper.limit)
    walkgraph$modularity = modularity(walkgraph)
    walkgraph$membership = membership(walkgraph)
    return(walkgraph)
  } else if (method=="iterative_spinglass")

      if (signed==TRUE){
        cores = detectCores()
        cl = makeCluster(cores-1)
        clusterExport(cl, c("graph", "upper.limit", "signed", "repetitions"),envir=environment())
        clusterEvalQ(cl, library(igraph))
        print("Applying Repeated Spinglass Community Estimation")
        spins = parLapply(cl, 1:repetitions, function(r) {cluster_spinglass(graph, weights=E(graph)$weight, implementation = "neg", spins=upper.limit, update.rule = "config")})
        stopCluster(cl)
        Memberships = sapply(1:length(spins), function(s) spins[[s]]$membership)

      } else if (signed==FALSE){
        cores = detectCores()
        cl = makeCluster(cores-1)
        clusterExport(cl, c("graph", "upper.limit", "signed", "repetitions"),envir=environment())
        clusterEvalQ(cl, library(igraph))
        spins = parLapply(cl, 1:repetitions, function(r) {cluster_spinglass(graph, weights=abs(E(graph)$weight), spins=upper.limit, update.rule = "config")})
        stopCluster(cl)
        Memberships = sapply(1:length(spins), function(s) spins[[s]]$membership)
      }
      Memberships = t(Memberships)
      data=Memberships
      return(data)
  }
