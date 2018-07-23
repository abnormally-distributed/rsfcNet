#' Visualize a connectivity matrix
#' @param graph an igraph object or matrix
#' @param type  "static" (the default) creates a standard image, "html" creates an interactive image with plotly
#' @param save.html Save the HTML file (only for type="html")
#' @param plot.title A title for the image
#' @param directory Where to save the HTML for a plotly image. Deaults to current working directory.
#' @param filename What to call the html file
#' @param color.scheme A color scheme. Defaults to viridis from the viridis package, a colorblind friendly palette.
#' @export
#'
#' @details This creates a visualization of a connectivity matrix. Either a standard static image
#' can be created or an interactive image with plotly can be created. Color palettes can be supplied
#' that will work with either output type with packages such as viridis, colorRamps, grDevices, and
#' RColorBrewer.
#'
#' A few color palettes I recommend besides the default:
#'
#' colorRamps::blue2red(500) , which is particuarly good for visualizing a contrast
#' between positive and negative correlations
#'
#' colorRamps::ygobb(500))
#'
#' colorRamps::matlab.like(500)
#'
#' colorRamps::matlab.like2(500)
#'
#' viridis::magma(500)
#'
#' viridis::plasma(500)
#'
#' grDevices::grey.colors(500)
#'
#' @author Brandon Vaughan
#'
#' This creates a correlogram (heatmap) with plotly in html format.
heatmap = function(graph, plot=TRUE, type="static" ,save.html = FALSE, plot.title="Network Heatmap",
        directory=paste0(getwd(),"/"), filename= "heatmap.html", color.scheme=viridis::viridis(500)) {
  if (is.igraph(graph) == "TRUE") {
    graph = as.matrix(igraph::as_adjacency_matrix(graph,
                                                  edges = FALSE, attr = "weight", sparse = TRUE))
  } else if (is.matrix(graph)=="TRUE"){
    graph=graph
  }

  if (type=="html") {

  vals <- unique(scales::rescale(c(graph)))
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric(color.scheme, domain = NULL)(vals)
  colz <- setNames(data.frame(vals[o], cols[o]), NULL)
  p <- plotly::plot_ly(z = graph, colorscale = colz, type = "heatmap")  %>%
    plotly::layout(
      title = plot.title)

  if (plot==TRUE){
    return(p)
  }
  if (save.html==TRUE) {
    htmlwidgets::saveWidget(as_widget(p), file = paste0(directory,filename), selfcontained=TRUE)
  }

  } else if (type=="static") {

    if(!is.null(dev.list())) dev.off()
    fields::image.plot(graph, col=color.scheme, axes=FALSE)
    title(plot.title)
  }

}
