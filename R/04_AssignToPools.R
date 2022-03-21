#' Create the graph pooling filter, given the adjacency matrix of the input graph.
#' @include clusteringfunctions.R
#' @param modelInput An object of type "ModelInput".
#' @param method One of "kmeans", "hierarchical", "supernodes", or "betweenness".
#' @param poolCount Number of pools. If more than one layer of pooling is desired,
#' this parameter should be a vector where the first number is the number of
#' pools in the first layer.
#' @param
#' @export
DoPooling <- function(modelInput, method = "kmeans", layers = 1, poolCount){
  
  # Check that the pool count makes sense.
  if(min(poolCount) <= 1){
    stop("Each layer must have more than one pool.")
  }
  else if(max(poolCount) > length(modelInput@true.values)){
    stop("Number of pools cannot be greater than number of samples.")
  }
  else{
    for(i in 1:length(poolCount)-1){
      if(poolCount[i] < poolCount[i+1]){
        stop("Number of pools must decrease in subsequent layers.")
      }
    }
    
    # Find first pooling layer.
    pools <- list()
    firstPoolCount <- poolCount
    if(length(poolCount) == 1){
      firstPoolCount <- poolCount[1]
    }
    if(method == "hierarchical"){
      pools[[1]] <- CreateHierarchicalPoolingFilter(modelInput = modelInput, 
                                                    poolCount = firstPoolCount)
    }
    if(method == "kmeans"){
      #pools[[1]] <- CreateKmeansPoolingFilter(modelInput = modelInput, 
      #                                              poolCount = firstPoolCount)
    }
  }
}

#' Create the graph pooling filter, given the adjacency matrix of the input graph.
#' @include clusteringfunctions.R
#' @param modelInputs An object of type "ModelInput".
#' @param poolCount Number of pools. 
#' @export
CreateHierarchicalPoolingFilter <- function(modelInputs, poolCount){
  
  # Perform hierarchical clustering.
  hier <- doHierarchicalClustering(modelInputs)
  
  # Arrange matrix of mappings.
  clusters <- stats::cutree(hier, k = poolCount)[rownames(modelInputs@node.wise.prediction)]
  cluster_count <- length(unique(clusters))
  mappings <- matrix(0, ncol = cluster_count, nrow = nrow(modelInputs@node.wise.prediction))
  cluster_length <- matrix(0, cluster_count)
  for(i in 1:cluster_count){
    which_in_cluster <- which(clusters == i)
    cluster_length[i] <- length(which_in_cluster)
    mappings[which_in_cluster, i] <- 1
  }
  newPoolingFilter <- methods::new("PoolingFilter", filter=mappings, filter.type="none",
                                   cluster.sizes=cluster_length, individual.filters=list())
  return(newPoolingFilter)
}