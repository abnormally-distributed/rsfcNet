#' Create Gephi File from graph object
#'
#' @param graph An igraph object
#' @param node_names A vector of node labels. By default searches global environemnt for a vector titled "colnames".
#' @param node_attributes A vector or data frame of different metrics describing the nodes (ie, modularity, degree).
#' @param directory Character vector specifyig directory for file output. Defaults to current working directory.
#' @param name A name for the output. Defaults to the current date.
#'
#' @details This function simplifies the use of the igraph to .gexf conversion in the rgexf package.
#'
#' @export
#' @author Brandon Vaughan
#'
#' @examples
#' attributes = cbind.data.frame(colMeans(EVC) , colMean(STR), colMeans(MODULE))
#' write_gephi(graph=group_avg_graph, attributes, directory = "C:/Users/Abnormally-Distributed/Documents/", name="Group_Average_Graph")
#'
#' @references
#'
#' Bastian M., Heymann S., Jacomy M. (2009). Gephi: an open source software for exploring and manipulating networks. International AAAI Conference on Weblogs and Social Media.
#'
#' https://gephi.org/
#
write_gephi = function(graph, node_attributes, directory=getwd(), name=Sys.Date(), node_names = NULL) {
  nodes_df = data.frame(ID = c(1:vcount(graph)), NAME = node_names)
  edges_att = data.frame(EDGE.WEIGHT= E(graph)$weight)
  edges_df = as.data.frame(get.edges(graph, c(1:ecount(graph))))
  nodes_coord_a = as.data.frame(layout.auto(graph, weights = E(graph)$weight, niter = 50))
  nodes_coord = cbind(nodes_coord_a, rep(0, times = nrow(nodes_coord_a)))
  file_title = paste0(directory,name,".gexf",sep="")
  print("Converting to gephi format. This may take 3-5 minutes depending on network size and computer speed.")
  write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = node_attributes, edgesAtt = edges_att, defaultedgetype = "undirected", output = file_title)
}
