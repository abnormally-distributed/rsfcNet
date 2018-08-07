#' Calculate the nodewise-deletion change in graph energy for a single graph
#'
#' This function calculates the change in graph energy observed due to the deletion of each node in the graph.
#' @param graph A network as an igraph object
#' @return A vector of centrality scores for each node.
#' @export
#' @author Brandon Vaughan
#'
#' @details Fornito, Zalesky, and Bullmore define a delta centrality measure as one that measures the change in
#' a global property of the graph that occurs due to the deletion of a node or edge (Fornito et al, 2016).
#' One informative global property is the energy of a graph. Graph energy was originally
#' applied in  organic chemistry to quantify the stability of molecular orbitals
#' associated with pi-electrons (Li, Shi, & Gutman, 2012). The graph energy informs about the connectivity of the graph
#' as a whole. A graph energy of zero means the nodes are not connected at all. Therefore,
#' to understand the importance of a particular node to the connectivity of the graph
#' a node can be removed and the change in graph energy used as a measure of node importance.
#'
#'
#' @seealso
#' \code{\link[rsfcNet]{laplace_centr}}
#' \code{\link[rsfcNet]{graph_energy}}
#' \code{\link[igraph]{fiedler_value}}
#'
#' @examples
#'
#' delta_e = delta_energy(Graph)
#'
#' @references
#'  Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Li, X.; Shi, Y.; Gutman, I. (2012), Graph Energy, New York: Springer, ISBN 978-1-4614-4219-6.

delta_energy = function(graph) {

  if (is.igraph(graph) == "TRUE") {
    graph = as.matrix(igraph::as_adjacency_matrix(graph, edges = FALSE, attr = "weight", sparse = TRUE))
  }

  diag(graph) <- 0
  graph_energy = sum(abs(eigen(graph, only.values=TRUE)$values))

  igraph = igraph::graph.adjacency(graph, mode="undirected", diag=FALSE, weighted=T)
  nodes = V(igraph)

  energy_after_drop <- c()
  for (i in V(igraph)) {
    g2 = delete_vertices(igraph, nodes[i])
    mat = as.matrix(igraph::as_adjacency_matrix(g2,  edges = FALSE, attr = "weight", sparse = TRUE))
    energy_after_drop[i] <- sum(abs(eigen(mat, only.values=TRUE)$values))
  }

  delta_energy = graph_energy-energy_after_drop
  return(delta_energy)
}

#' Calculate the nodewise-deletion change in graph energy for a list of graphs
#'
#' This function calculates the change in graph energy observed due to the deletion of each node in the graph.
#' @param graphs A list of igraph objects
#' @return A matrix of centrality scores for each node in each graph.
#' @export
#' @author Brandon Vaughan
#'
#' @details Fornito, Zalesky, and Bullmore define a delta centrality measure as one that measures the change in
#' a global property of the graph that occurs due to the deletion of a node or edge (Fornito et al, 2016).
#' One informative global property is the energy of a graph. Graph energy was originally
#' applied in  organic chemistry to quantify the stability of molecular orbitals
#' associated with pi-electrons (Li, Shi, & Gutman, 2012). The graph energy informs about the connectivity of the graph
#' as a whole. A graph energy of zero means the nodes are not connected at all. Therefore,
#' to understand the importance of a particular node to the connectivity of the graph
#' a node can be removed and the change in graph energy used as a measure of node importance.
#'
#'
#' @seealso
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#' \code{\link[rsfcNet]{graph_energy_mult}}
#' \code{\link[igraph]{fiedler_value_mult}}
#'
#' @examples
#'
#' delta_e = delta_energy_mult(Graphs)
#'
#' @references
#'  Fornito, A., Zalesky, A., & Bullmore, E. (2016). Centrality and Hubs. Chapter 5. Fundamentals of Brain Network Analysis, 137-161. doi:10.1016/b978-0-12-407908-3.00005-4
#'
#' Li, X.; Shi, Y.; Gutman, I. (2012), Graph Energy, New York: Springer, ISBN 978-1-4614-4219-6.


delta_energy_mult = function(graphs) {
  sapply(graphs, rsfcNet::delta_energy)
}

#' Calculate the graph energy for a single graph
#'
#' This function calculates the energy of a graph.
#' @param graph A network as an igraph object or connectivity matrix
#' @return A numeric value of the graph energy
#' @export
#' @author Brandon Vaughan
#'
#' @details Graph energy was originally
#' applied in  organic chemistry to quantify the stability of molecular orbitals
#' associated with pi-electrons (Li, Shi, & Gutman, 2012). The graph energy informs about the connectivity of the graph
#' as a whole. A graph energy of zero means the nodes are not connected at all. The graph energy
#' is calculated by summing the absolute values of the eigenvalues for the matrix. See Daianu et al
#' (2015) for an example of an application to neuroimaging.
#'
#'
#' @seealso
#' \code{\link[rsfcNet]{laplace_centr}}
#' \code{\link[rsfcNet]{delta_energy}}
#' \code{\link[igraph]{fiedler_value}}
#' \code{\link[igraph]{graph_energy_mult}}
#'
#' @examples
#'
#' energy = graph_energy(Graph)
#'
#' @references
#' Daianu, M., Mezher, A., Jahanshad, N., Hibar, D. P., Nir, T. M., Jack, C. R., … Thompson, P. M. (2015). SPECTRAL GRAPH THEORY AND GRAPH ENERGY METRICS SHOW EVIDENCE FOR THE ALZHEIMER’S DISEASE DISCONNECTION SYNDROME IN APOE-4 RISK GENE CARRIERS. IEEE International Symposium on Biomedical Imaging, 2015, 458–461. http://doi.org/10.1109/ISBI.2015.7163910
#'
#' Li, X.; Shi, Y.; Gutman, I. (2012), Graph Energy, New York: Springer, ISBN 978-1-4614-4219-6.

graph_energy = function(graph) {
  if (is.igraph(graph)=="TRUE"){
    graph = as.matrix(as_adj(graph, type="both", edges=FALSE, sparse=TRUE))
  }
  eigen = eigen(graph)$values
  sum(abs(eigen))
}

#' Calculate the graph energy for a list of graphs
#'
#' This function calculates the energy each graph in a list.
#' @param graphs A list of networks as igraph objects or as matrices
#' @return A vector of graph energy values for each graph.
#' @export
#' @author Brandon Vaughan
#'
#' @details Graph energy was originally
#' applied in  organic chemistry to quantify the stability of molecular orbitals
#' associated with pi-electrons (Li, Shi, & Gutman, 2012). The graph energy informs about the connectivity of the graph
#' as a whole. A graph energy of zero means the nodes are not connected at all. The graph energy
#' is calculated by summing the absolute values of the eigenvalues for the matrix. See Daianu et al
#' (2015) for an example of an application to neuroimaging.
#'
#'
#' @seealso
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#' \code{\link[rsfcNet]{delta_energy_mult}}
#' \code{\link[igraph]{fiedler_value_mult}}
#'
#' @examples
#'
#' energy = graph_energy_mult(Graphs)
#'
#' @references
#' Daianu, M., Mezher, A., Jahanshad, N., Hibar, D. P., Nir, T. M., Jack, C. R., … Thompson, P. M. (2015). SPECTRAL GRAPH THEORY AND GRAPH ENERGY METRICS SHOW EVIDENCE FOR THE ALZHEIMER’S DISEASE DISCONNECTION SYNDROME IN APOE-4 RISK GENE CARRIERS. IEEE International Symposium on Biomedical Imaging, 2015, 458–461. http://doi.org/10.1109/ISBI.2015.7163910
#'
#' Li, X.; Shi, Y.; Gutman, I. (2012), Graph Energy, New York: Springer, ISBN 978-1-4614-4219-6.

graph_energy_mult = function(graphs) {

  sapply(graphs, rsfcNet::graph_energy)

}


#' Calculate the fiedler value for a graph
#'
#' This function calculates the fiedler value for a graph.
#' @param graph A network as an igraph object
#' @return A numeric value of the fiedler value.
#' @export
#' @author Brandon Vaughan
#'
#' @details The Fiedler value is the second smallest eigenvalue of the laplacian representation
#' of a graph. The closer the Fiedler value is to zero the more easily the graph can be split into
#' separate components unconnected to each other. The Fiedler value is also known as the algebraic
#' connectivity of a graph (Mohar, 1991). Hence the fiedler value can be used as a measure of a network's
#' robustness to becoming disconnected. See Daianu et al (2014) for an application to neuroimaging.
#'
#'
#' @seealso
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#' \code{\link[rsfcNet]{delta_energy_mult}}
#' \code{\link[igraph]{fiedler_value_mult}}
#'
#' @examples
#'
#' energy = fiedler_value(Graph)
#'
#' @references
#' Mohar, B., The Laplacian Spectrum of Graphs, in Graph Theory, Combinatorics, and Applications, Vol. 2, Ed. Y. Alavi, G. Chartrand, O. R. Oellermann, A. J. Schwenk, Wiley, 1991, pp. 871–898.
#'
#' Daianu, M., Jahanshad, N., Nir, T. M., Leonardo, C. D., Jack, C. R., Weiner, M. W., … Thompson, P. M. (2014). Algebraic connectivity of brain networks shows patterns of segregation leading to reduced network robustness in Alzheimer’s disease. Computational Diffusion MRI: MICCAI Workshop, Boston, MA, USA, September 2014.  MICCAI Workshop on Computation Diffusion MRI (Boston, MA), 55–64. http://doi.org/10.1007/978-3-319-11182-7_6

fiedler_value = function(graph) {

  graph = abs(as.matrix(as_adj(graph, attr="weight", type="both", edges=FALSE)))
  graph = graph.adjacency(graph, weighted=T, mode="undirected", diag=FALSE)
  laplacian = as.matrix(graph.laplacian(graph))
  eigen = eigen(laplacian)$values
  n <- length(eigen)
  eigen[n-1]

}


#' Calculate the fiedler value for a list of graphs.
#'
#' This function calculates the fiedler value for a graph.
#' @param graphs A list of igraph objects
#' @return A numeric value of the fiedler value.
#' @export
#' @author Brandon Vaughan
#'
#' @details The Fiedler value is the second smallest eigenvalue of the laplacian representation
#' of a graph. The closer the Fiedler value is to zero the more easily the graph can be split into
#' separate components unconnected to each other. The Fiedler value is also known as the algebraic
#' connectivity of a graph (Mohar, 1991). Hence the fiedler value can be used as a measure of a network's
#' robustness to becoming disconnected. See Daianu et al (2014) for an application to neuroimaging.
#'
#'
#' @seealso
#' \code{\link[rsfcNet]{laplace_centr_mult}}
#' \code{\link[rsfcNet]{delta_energy_mult}}
#' \code{\link[igraph]{graph_energy_mult}}
#' \code{\link[igraph]{fiedler_value}}
#'
#' @examples
#'
#' energy = fiedler_value_mult(Graphs)
#'
#' @references
#' Mohar, B., The Laplacian Spectrum of Graphs, in Graph Theory, Combinatorics, and Applications, Vol. 2, Ed. Y. Alavi, G. Chartrand, O. R. Oellermann, A. J. Schwenk, Wiley, 1991, pp. 871–898.
#'
#' Daianu, M., Jahanshad, N., Nir, T. M., Leonardo, C. D., Jack, C. R., Weiner, M. W., … Thompson, P. M. (2014). Algebraic connectivity of brain networks shows patterns of segregation leading to reduced network robustness in Alzheimer’s disease. Computational Diffusion MRI: MICCAI Workshop, Boston, MA, USA, September 2014.  MICCAI Workshop on Computation Diffusion MRI (Boston, MA), 55–64. http://doi.org/10.1007/978-3-319-11182-7_6

fiedler_value_mult = function(graphs) {

sapply(graphs, rsfcNet::fiedler_value)

}

