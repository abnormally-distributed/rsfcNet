#' Reiterative louvain community detection for signed networks
#'
#' Get the consensus community memberships from multiple iterations of a louvain algorithm
#'
#' @param graph an igraph object
#' @param iter the number of iterations. defaults to 100 but more may be desired.
#' @export
#' @author Brandon Vaughan
#' @details Gets the consensus for modules
#' @examples
#' **##Not run**
#' iterative_louvain_signed(graph)
#'  ## End(**Not run**)
#'
iterative_louvain_signed <-function(graph,iter=100,cores=1)
{
  original.graph = graph
  original.iter = iter
  graph = as.matrix(igraph::as_adjacency_matrix(graph,
                                                edges = FALSE, attr = "weight", sparse = TRUE))
  graph = graph * (graph>0)
  graph = igraph::graph.adjacency(graph, weighted = T, mode = "undirected", diag = TRUE)
  nodes <- as.vector(seq(1:length(V(graph))))
  edgelist <- visNetwork::toVisNetworkData(graph)$edge
  edgelist$weight <- abs(edgelist$weight)
  iter = round(iter * .75)

  if(is.data.frame(edgelist) & ncol(edgelist)==3)
  {
    colnames(edgelist)<-c("V1","V2","weight")
    row.names(edgelist)<-seq(1,nrow(edgelist))
    if(cores==1)
    {
      results<-foreach(1:iter,.combine=cbind) %do% {
        rand_edge<-edgelist[order(rnorm(nrow(edgelist))),]
        Bn<-igraph::graph_from_data_frame(rand_edge,directed = F)
        brain_community<-igraph::cluster_louvain(Bn)
        mod_members<-igraph::membership(brain_community)
        mod_members[order(as.numeric(names(mod_members)))]

      }
    } else if (cores >1 )
    {
      cl<-makeCluster(cores)
      registerDoParallel(cl)
      results<-foreach(1:iter,.combine=cbind) %dopar% {
        rand_edge<-edgelist[order(rnorm(nrow(edgelist))),]
        Bn<-igraph::graph_from_data_frame(rand_edge,directed = F)
        brain_community<-igraph::cluster_louvain(Bn)
        mod_members<-igraph::membership(brain_community)
        mod_members[order(as.numeric(names(mod_members)))]

      }
      stopCluster(cl)
    } else {
      stop("doParallel and foreach not working")
    }
    merged_results.pos<-merge(as.data.frame(nodes),results,by.x="nodes",by.y="row.names",all.x=TRUE)
    merged_results.pos[is.na(merged_results.pos)]<-0
  }

  Qpos = rsfcNet:::consensus.utility.func(merged_results.pos)
  Qpos = modularity(graph, Qpos$membership, weights = abs(E(graph)$weight))
  v.pos = sum(strength_signed(graph)$positive_strength)

  ################################# Negative Modularity ##############################################

  graph = original.graph
  graph = as.matrix(igraph::as_adjacency_matrix(graph,
                                                edges = FALSE, attr = "weight", sparse = TRUE))
  graph = graph * (graph < 0)
  graph = igraph::graph.adjacency(graph, weighted = T, mode = "undirected", diag = TRUE)
  nodes <- as.vector(seq(1:length(V(graph))))
  edgelist <- visNetwork::toVisNetworkData(graph)$edge
  edgelist$weight <- abs(edgelist$weight)
  iter = round(original.iter * .25)

  if(is.data.frame(edgelist) & ncol(edgelist)==3)
  {
    colnames(edgelist)<-c("V1","V2","weight")
    row.names(edgelist)<-seq(1,nrow(edgelist))
    if(cores==1)
    {
      results<-foreach(1:iter,.combine=cbind) %do% {
        rand_edge<-edgelist[order(rnorm(nrow(edgelist))),]
        Bn<-igraph::graph_from_data_frame(rand_edge,directed = F)
        brain_community<-igraph::cluster_louvain(Bn)
        mod_members<-igraph::membership(brain_community)
        mod_members[order(as.numeric(names(mod_members)))]

      }
    } else if (cores >1 )
    {
      cl<-makeCluster(cores)
      registerDoParallel(cl)
      results<-foreach(1:iter,.combine=cbind) %dopar% {
        rand_edge<-edgelist[order(rnorm(nrow(edgelist))),]
        Bn<-igraph::graph_from_data_frame(rand_edge,directed = F)
        brain_community<-igraph::cluster_louvain(Bn)
        mod_members<-igraph::membership(brain_community)
        mod_members[order(as.numeric(names(mod_members)))]

      }
      stopCluster(cl)
    } else {
      stop("doParallel and foreach not working")
    }
    merged_results.neg<-merge(as.data.frame(nodes),results,by.x="nodes",by.y="row.names",all.x=TRUE)
    merged_results.neg[is.na(merged_results.neg)]<-0

    Qneg= rsfcNet:::consensus.utility.func(merged_results.neg)
    Qneg = -modularity(graph, Qneg$membership, weights = abs(E(graph)$weight))
    v.neg = -1 * sum(strength_signed(graph)$positive_strength) ## Because the negative graph is absolute valued

    Qstar = Qpos + ((v.neg / (v.neg+v.pos))*Qneg)
    Qsimple = Qpos + Qneg

    #################### Combine partitions ###########################################

    merged_results = cbind(merged_results.pos, merged_results.neg)

    modules = t(merged_results[,-1])
    node.count = ncol(modules)
    node.names <- as.character(seq(1:node.count))
    colnames(modules) = node.names
    output = list()
    for (i in 1:nrow(modules)) {
      row = as.numeric(modules[i, ])
      output.inner = list()
      for (j in 1:length(row)) {
        module = character()
        module = c(module, colnames(modules)[which(row ==
                                                     row[j])])
        output.inner[[j]] = module
      }
      output.inner = unique(output.inner)
      output[[i]] = output.inner
    }

    consensus_num_modules = as.numeric(names(sort(table(unlist(lapply(output,
                                                                      length))), decreasing = TRUE))[1])
    output = output[lapply(output, length) == consensus_num_modules]
    module = list()
    for (i in 1:consensus_num_modules) {
      list.common_elements = list()
      for (p in 1:length(output)) {
        list.common_elements[[p]] = unlist(output[[p]][i])
      }
      common_elements = Reduce(intersect, list.common_elements)
      module[[i]] = common_elements
      for (p in 1:length(list.common_elements)) {
        vector = setdiff(list.common_elements[[p]], common_elements)
        if (length(vector) > 0) {
          for (j in 1:length(vector)) {
            counter = vector(length = length(list.common_elements))
            for (k in 1:length(list.common_elements)) {
              counter[k] = vector[j] %in% list.common_elements[[k]]
            }
            if (length(which(counter == TRUE)) >= ceiling((length(counter)/2) +
                                                          0.001)) {
              module[[i]] = c(module[[i]], vector[j])
            }
          }
        }
      }
    }
    module = lapply(module, unique)
    wanderers = setdiff(colnames(modules), unlist(module))
    if (length(wanderers) > 0) {
      for (pp in 1:length(wanderers)) {
        temp = matrix(nrow = length(output), ncol = consensus_num_modules)
        for (i in 1:nrow(temp)) {
          for (j in 1:ncol(temp)) {
            temp[i, j] = wanderers[pp] %in% unlist(output[[i]][j])
          }
        }
        index.best = which(temp[1, ] == TRUE)
        module[[index.best]] = c(module[[index.best]], wanderers[pp])
      }
    }
    module = lapply(module, function(m) as.numeric(m))
    membership = module
    membership = unlist(lapply(1:length(membership), function(i) rep(i,
                                                                     length(membership[[i]]))))
    membership = cbind.data.frame(unlist(module), membership)
    membership = membership$membership[order(membership$`unlist(module)`)]
    Module = list(groups = module, membership = membership)
    module = list(groups = Module$groups, membership=Module$membership, Qsimple=Qsimple, Qstar=Qstar)
    return(module)
  }
  else {
    stop("Edges are not a data frame or does not have 3 columns")
  }
}


#' Utility function for louvain signed algorithm
#'
#' Get the consensus community memberships from multiple iterations of a louvain algorithm
#'
#' @param merged_results the merged results
#' @param iter the number of iterations. defaults to 100 but more may be desired.
#' @keywords internal
#' @author Brandon Vaughan
#' @details Gets the consensus for modules
#' @examples
#' **##Not run**
#' consensus.utility.func(merged_results)
#'  ## End(**Not run**)
#'
#'
#'
consensus.utility.func = function(merged_results) {

  modules = t(merged_results[,-1])
  node.count = ncol(modules)
  node.names <- as.character(seq(1:node.count))
  colnames(modules) = node.names
  output = list()
  for (i in 1:nrow(modules)) {
    row = as.numeric(modules[i, ])
    output.inner = list()
    for (j in 1:length(row)) {
      module = character()
      module = c(module, colnames(modules)[which(row ==
                                                   row[j])])
      output.inner[[j]] = module
    }
    output.inner = unique(output.inner)
    output[[i]] = output.inner
  }

  consensus_num_modules = as.numeric(names(sort(table(unlist(lapply(output,
                                                                    length))), decreasing = TRUE))[1])
  output = output[lapply(output, length) == consensus_num_modules]
  module = list()
  for (i in 1:consensus_num_modules) {
    list.common_elements = list()
    for (p in 1:length(output)) {
      list.common_elements[[p]] = unlist(output[[p]][i])
    }
    common_elements = Reduce(intersect, list.common_elements)
    module[[i]] = common_elements
    for (p in 1:length(list.common_elements)) {
      vector = setdiff(list.common_elements[[p]], common_elements)
      if (length(vector) > 0) {
        for (j in 1:length(vector)) {
          counter = vector(length = length(list.common_elements))
          for (k in 1:length(list.common_elements)) {
            counter[k] = vector[j] %in% list.common_elements[[k]]
          }
          if (length(which(counter == TRUE)) >= ceiling((length(counter)/2) +
                                                        0.001)) {
            module[[i]] = c(module[[i]], vector[j])
          }
        }
      }
    }
  }
  module = lapply(module, unique)
  wanderers = setdiff(colnames(modules), unlist(module))
  if (length(wanderers) > 0) {
    for (pp in 1:length(wanderers)) {
      temp = matrix(nrow = length(output), ncol = consensus_num_modules)
      for (i in 1:nrow(temp)) {
        for (j in 1:ncol(temp)) {
          temp[i, j] = wanderers[pp] %in% unlist(output[[i]][j])
        }
      }
      index.best = which(temp[1, ] == TRUE)
      module[[index.best]] = c(module[[index.best]], wanderers[pp])
    }
  }
  module = lapply(module, function(m) as.numeric(m))
  membership = module
  membership = unlist(lapply(1:length(membership), function(i) rep(i,
                                                                   length(membership[[i]]))))
  membership = cbind.data.frame(unlist(module), membership)
  membership = membership$membership[order(membership$`unlist(module)`)]
  module = list(groups = module, membership = membership)
  return(module)
}


