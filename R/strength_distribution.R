#' Get the degree or strength distribution for a graph
#'
#' Get either the degree distribution (as a count or as a density) or the strength distribution for a graph.
#'
#' @param graph A network as an igraph object or connectivity matrix.
#' @param cumulative Should the cumulative sum be returned or the density? Defaults to cumulative=FALSE
#' @param signed If TRUE (the default) uses the strength_signed function to calculate the distribution.
#' @param which.sign If "star" returns the weighted strength* metric, if "pos" returns the distribution of the
#' positive connections, and if "neg" returns the distribution of the negative connections.
#' @param count If the weights are binary can be set to count=TRUE to get the counts of each degree instead
#' of the probability densities.
#'
#' @author Brandon Vaughan
#' @export
#'
#' @examples
#' #' **## Not run:**
#' # If you want to examine the distributions of positive and negative connections
#' pos = strength_distribution(graph, which.sign="pos")
#' neg = strength_distribution(graph, which.sign="neg")
#'
#' ## End(**Not run**)
#'
#' @seealso
#' \code{\link[rsfcNet]{dist_fit}}
#' \code{\link[rsfcNet]{strength_signed}}
#'
strength_distribution <- function (graph, cumulative = FALSE, signed=TRUE, which.sign="star", count=FALSE)
{
  if (igraph::is.igraph(graph) == FALSE) {
    graph = igraph::graph.adjacency(graph, mode = "undirected",
                                    weighted = T, diag = FALSE)
  }

  testInteger <- function(x){
    test <- all(x == as.integer(x))
    if(test == TRUE){ return(TRUE) }
    else { return(FALSE) }
  }

  if(testInteger(as.vector(as_adj(graph, attr="weight")))=="TRUE") {
    cs <- degree(graph)
    if (count=="TRUE"){
      hi <- hist(cs, -1:max(cs), plot = FALSE)$count
    } else {
      hi <- hist(cs, -1:max(cs), plot = FALSE)$density}
    if (!cumulative) {
      res <- hi
    }
    else {
      res <- rev(cumsum(rev(hi)))
    }
    res
    return(res)
  } else {

    if (signed==FALSE) {
      cs <- graph.strength(graph, ...)
    } else if (signed==TRUE){
      if (which.sign=="pos") {
        cs <- strength_signed(graph, scale=FALSE)$positive_strength
      } else if (which.sign=="neg") {
        cs <- strength_signed(graph, scale=FALSE)$negative_strength
      } else if (which.sign=="star") {
        cs <- strength_signed(graph, scale=FALSE)$strength_star
      }
    }


    #value.range <- seq(from = min(cs)-1, to = max(cs)+1)
    hi <- density(cs)$y
    if (!cumulative) {
      res <- hi
    }
    else {
      res <- rev(cumsum(rev(hi)))
    }
    res
  }
}

