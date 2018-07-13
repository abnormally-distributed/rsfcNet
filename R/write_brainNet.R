#' Create Node and Edge Files for Brain Net Viewer
#'
#' @param MNIcoordinates A matrix with columns consisting of the X,Y,Z MNI152 coordinates (must be in this order) for ROIs
#' @param graph A network in either matrix format or as an igraph object
#' @param node.colors A vector of labels, such as module membership. Defaults to all the same color.
#' @param node.size A vector of numbers describing the nodes (such as degree)
#' @param directory Directory for file output. Defaults to current working directory.
#' @param name A name for the output. Defaults to the current date.
#'
#'
#'
#' @export
#' @author Brandon Vaughan
#'
#' @examples
#' brainnet.matrix = threshold_matrix(avg_matrix) # Threshold a weighted correlation matrix first
#' brainnet.matrix = threshold_and_binarize(avg_matrix) # Or binarize
#' write.brainNet(MNIcoordinates=Node_Atlas[,1:3], graph = brainnet.matrix, node.colors=module$membership, node.size=LEV.AVG, directory="C:/Users/Abnormally-Distributed/Documents/R/", name="leverage_network")
#'
#' @references
#'
#' Xia M, Wang J, He Y (2013) BrainNet Viewer: A Network Visualization Tool for Human Brain Connectomics. PLoS ONE 8: e68910.
#'
#'
#' https://www.nitrc.org/projects/bnv/
#'
write_brainNet = function(MNIcoordinates, graph, node.colors=rep(1, .GlobalEnv$n.nodes), node.size, node.labels=seq(1:.GlobalEnv$n.nodes), directory=getwd(), name=Sys.Date()){

  if (is.igraph(graph)=="TRUE"){
    graph = as.matrix(igraph::as_adjacency_matrix(graph, edges = FALSE, attr = "weight", sparse=TRUE)) }


  node.file = cbind.data.frame(MNIcoordinates, node.colors, node.size, node.labels)
  file.name.node = paste0(directory,name,".node",sep="")
  file.name.edge = paste0(directory,name,".edge",sep="")
  write.table(node.file, file=file.name.node, sep="\t", col.names = F, row.names = F)
  write.table(graph, file=file.name.edge, sep="\t", col.names = F, row.names = F)
}
