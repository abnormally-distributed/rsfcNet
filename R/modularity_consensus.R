#' Get the consensus community membership labels
#'
#' Get the consensus community memberships for multiple subjects or multiple iterations of an algorithm over a single graph.
#'
#' @param modules a matrix containing membership labels across each row for each node (represented by column).
#' @param node.names custom node names if desired
#' @export
#' @author Brandon Vaughan
#'
#' @details Gets a consensus on module membership for modules that have clear membership. The rest are left unsorted as "wanderers".
#' An overall crude method of getting consensus, but relaxes the assumption that all nodes in a network must belong to a module.
#'
#' @examples
#' **##Not run**
#' modularity_consensus(module_list, node.names=brain_regions)
#'  ## End(**Not run**)
#'
modularity_consensus =  function(modules, node.names="NA"){

  if (node.names=="NA"){
    node.count = ncol(modules)
    node.names <- as.character(seq(1:node.count))
  }

    output=list()
    for (i in 1:nrow(modules)){
      row=as.numeric(modules[i,])
      output.inner=list()
      for (j in 1:length(row)){
        module=character()
        module=c(module,colnames(modules)[which(row==row[j])])
        output.inner[[j]]=module
      }
      output.inner=unique(output.inner)
      output[[i]]=output.inner
    }

    # Get the  mode of the vector representing the number of modules found by each method
    consensus_num_modules=as.numeric(names(sort(table(unlist(lapply(output,length))),decreasing=TRUE))[1])

    # Removes nodes that do not correspond to this consensus
    output=output[lapply(output,length)==consensus_num_modules]

    # Now find nodes in a membership vector that belong to a common module and then get the majority
    # vote for nodes that are not part of the set

    module=list()

    for (i in 1:consensus_num_modules){
      list.common_elements=list()
      for (p in 1:length(output)){
        list.common_elements[[p]]=unlist(output[[p]][i])
      }

      # Candidate module_i
      common_elements=Reduce(intersect,list.common_elements)
      module[[i]]=common_elements

      # Reinforce the module.
      for (p in 1:length(list.common_elements)){
        vector=setdiff(list.common_elements[[p]],common_elements)
        if (length(vector)>0){
          for (j in 1:length(vector)){
            counter=vector(length=length(list.common_elements))
            for (k in 1:length(list.common_elements)){
              counter[k]=vector[j]%in%list.common_elements[[k]]
            }
            if(length(which(counter==TRUE))>=ceiling((length(counter)/2)+0.001)){
              module[[i]]=c(module[[i]],vector[j])
            }
          }
        }
      }
    }

    module=lapply(module,unique)

    # variables for which consensus has not been reached
    wanderers=setdiff(colnames(modules),unlist(module))

    if (length(wanderers)>0){
      for (pp  in 1:length(wanderers)){
        temp=matrix(nrow=length(output),ncol=consensus_num_modules)
        for (i in 1:nrow(temp)){
          for (j in 1:ncol(temp)){
            temp[i,j]=wanderers[pp]%in%unlist(output[[i]][j])
          }
        }
        # use the partition of the first method when no majority exists
        # (this allows ordering of partitions by decreasing modularity values for instance)
        index.best=which(temp[1,]==TRUE)
        module[[index.best]]=c(module[[index.best]],wanderers[pp])
      }
    }
    output=list(module=module,wanderers=wanderers)
    return(output)
}
