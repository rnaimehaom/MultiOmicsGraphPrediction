#' Create a line graph given the original, unweighted graph. In a line graph,
#' edges are nodes, and edges connected by a node are edges.
#' @param predictionsByEdge Prediction levels corresponding to each edge.
#' @param graphWithPredictions Original graph data frame.
#' @param edgeTypeList List containing one or more of the following:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
CreateLineGraph <- function(predictionsByEdge, graphWithPredictions, edgeTypeList){
  # Step 1: Define vertices.
  line_graph_vertices <- predictionsByEdge$Node
  list_to_add <- list()
  
  # Step 2: Identify edges that share "to" nodes.
  if("shared.outcome.analyte" %in% edgeTypeList){
    list_to_add[[length(list_to_add)+1]] <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                                                  graphWithPredictions = graphWithPredictions,
                                                                  nodeType1 = "to", nodeType2 = "to")
  }
  
  # Step 3: Identify edges that share "from" nodes.
  if("shared.independent.analyte" %in% edgeTypeList){
    list_to_add[[length(list_to_add)+1]] <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                                                  graphWithPredictions = graphWithPredictions,
                                                                  nodeType1 = "from", nodeType2 = "from")
  }
  
  # Step 4: For each "to" node, identify edges that share a "from" node.
  if("analyte.chain" %in% edgeTypeList){
    list_to_add[[length(list_to_add)+1]] <- FindEdgesSharingNodes(predictionsByEdge = predictionsByEdge,
                                                                  graphWithPredictions = graphWithPredictions,
                                                                  nodeType1 = "to", nodeType2 = "from")
  }
  
  # Step 6: Concatenate.
  shared_df <- list_to_add[[1]]
  if(length(list_to_add) > 1){
    shared_df <- do.call(rbind, list_to_add)
  }
  
  # Step 7: Convert to graph.
  line_graph <- igraph::graph_from_data_frame(shared_df)
  
  # Step 8: Convert to adjacency matrix.
  line_graph_A <- igraph::as_adj(line_graph)
  
  # Return.
  return(line_graph_A)
}

#' Find edges that share nodes and add them to a data frame.
#' @param predictionsByEdge Prediction levels corresponding to each edge.
#' @param graphWithPredictions Original graph data frame.
#' @param nodeType1 Either "to" or "from".
#' @param nodeType2 Either "to or "from".
FindEdgesSharingNodes <- function(predictionsByEdge, graphWithPredictions, nodeType1,
                                  nodeType2){
  # Convert predictions to data frame.
  graph_df <- igraph::as_data_frame(graphWithPredictions)
  
  # Find shared nodes.
  nodes <- unique(graph_df[,nodeType1])
  to_shared <- lapply(1:length(nodes), function(i){
    
    # Find all line graph nodes starting with or ending with the analyte in question.
    node <- nodes[i]
    combs <- NULL
    set_with_1 <- predictionsByEdge$Node[which(graph_df[,nodeType1] == node)]
    set_with_2 <- predictionsByEdge$Node[which(graph_df[,nodeType2] == node)]
    
    # If there are multiple line graph nodes including this analyte, return them.
    combs <- expand.grid(set_with_1, set_with_2)
    combs$Var1 <- as.character(combs$Var1)
    combs$Var2 <- as.character(combs$Var2)
    combs <- combs[which(combs$Var1 != combs$Var2),]
    colnames(combs) <- c("to", "from")
    return(combs)
  })
  
  # Concatenate all shared nodes.
  line_graph_df <- do.call(rbind, to_shared)
  
  return(line_graph_df)
}

#' Predict Y given current weights.
#' @param modelResults An object of the ModelResults class.
#' @param iteration The current iteration.
#' @param pooling Whether or not to pool the weights.
#' @param convolution Whether or not to perform convolution.
DoPrediction <- function(modelResults, iteration, pooling, convolution){
  # Propagate forward.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  Y.pred <- X
  # If convolution is to be performed, perform convolution.
  if(convolution == TRUE){
    Y.pred.list <- lapply(1:dim(X)[2], function(i){
      return(A.hat %*% X[,i])
    })
    Y.pred <- do.call(cbind, Y.pred.list)
  }
  # If pooling is to be performed, multiply by pool either after or before weights
  # are learned, respectively.
  if(pooling == TRUE){
    if(modelResults@weights.after.pooling == TRUE){
      S_all <- modelResults@pooling.filter@individual.filters
      if(modelResults@current.iteration == 1){
        S_all <- CreateFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                              X = X)
        modelResults@pooling.filter@individual.filters <- S_all
      }
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(Y.pred[,i]) %*% S_all[[i]] * Theta.old[,i]))
      }))
    }else{
      S_all <- AdjustFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                            X = X, Theta = Theta.old)
      modelResults@pooling.filter@individual.filters <- S_all
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(Y.pred[,i] * Theta.old[,i]) %*% S_all[[i]]))
      }))
    }
  }else{
    Y.pred <- unlist(lapply(1:(dim(Theta.old)[2]), function(i){
      return(sum(t(Y.pred[,i] * Theta.old[,i])))
    }))
  }
  
  # Use activation function if output is of a character type. Note that all character
  # types are converted into factors, and since only binary factors are accepted by
  # the package, the values will be 1 (for the alphanumerically lowest level) and 2
  # (for the alphanumerically highest level).
  if(modelResults@model.input@outcome.type == "categorical"){
    if(modelResults@activation.type == "softmax"){
      Y.pred <- round(SoftmaxWithCorrection(Y.pred))
    }else if(modelResults@activation.type == "tanh"){
      Y.pred <- round(TanhWithCorrection(Y.pred))
    }else{
      Y.pred <- round(SigmoidWithCorrection(Y.pred))
    }
  }else{
    Y.pred <- Y.pred / dim(Theta.old)[1]
  }
  
  # Return the prediction.
  return(Y.pred)
}

#' Train the graph learning model, using the specifications in the ModelResults
#' class and storing the results in the ModelResults class.T
#' @param modelResults An object of the ModelResults class.
#' @param iteration The current iteration.
#' @param pooling Whether or not to pool the weights.
#' @param convolution Whether or not to perform convolution.
#' @param ridgeRegressionWeight The hyperparameter weight assigned
#' to the ridge regression parameter (often referred to as lambda in the
#' literature)
#' @param varianceWeight The hyperparameter weight assigned to the difference
#' in variances between Y and the predicted value of Y.
DoSingleTrainingIteration <- function(modelResults, iteration, pooling, convolution,
                                      ridgeRegressionWeight, varianceWeight){
  # Predict Y.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.weights, dim(X)[2]), ncol = dim(X)[2])
  Y.pred <- DoPrediction(modelResults, iteration, pooling, convolution)
  modelResults@outcome.prediction <- Y.pred
  
  # Backpropagate and calculate the error.
  Theta.new <- Theta.old
  error <- modelResults@iteration.tracking$Error[iteration-1]
  modelResults <- BackpropagateSingleLayer(modelResults, iteration, convolution,
                                           pooling, ridgeRegressionWeight,
                                           varianceWeight)
  if(modelResults@model.input@outcome.type == "categorical"){
    modelResults@iteration.tracking$Error[iteration+1] <- 
      ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
  }else{
    modelResults@iteration.tracking$Error[iteration+1] <- 
      ComputeNRMSE(modelResults@model.input@true.phenotypes, Y.pred)
  }
  
  # Modify the model results and return.
  return(modelResults)
}

#' Adjust the pooling filter (necessary for median, min, or max pooling).
#' @param poolingFilter The original pooling filter for the model, contained
#' within the modelResults object.
#' @param A.hat The normalized Laplacian for convolution.
#' @param X The input predictions.
AdjustFilter <- function(poolingFilter, A.hat, X){
  S <- poolingFilter@filter
  S_all <- lapply(1:dim(X)[2], function(i){
    return(S)
  })
  
  # For mean filters, divide each column by the count of features and
  # replicate for each sample.
  if(poolingFilter@filter.type == "mean"){
    S_all <- lapply(1:dim(X)[2], function(i){
      return(S / colSums(S))
    })
    # For min, max, or median filters, select the feature meeting the criteria
    # for each sample and return.
  }else if(poolingFilter@filter.type == "median"){
    S_all <- lapply(1:dim(X)[2], function(i){
      to_mult <- t(A.hat %*% X[,i])
      S.copy <- matrix(0, nrow = dim(S)[1], ncol = dim(S)[2])
      for(i in 1:dim(S)[2]){
        to_mult_filt <- to_mult[which(S[,i]!=0)]
        
        # Get median.
        if (length(to_mult_filt) %% 2 != 0) {
          sel_samp <- which(S[,i]!=0)[which(to_mult_filt == stats::median(to_mult_filt))[1]]
        } else if (length(to_mult_filt) %% 2 == 0) {
          a = sort(to_mult_filt)[c(length(to_mult_filt)/2, length(to_mult_filt)/2+1)]
          sel_samp <- which(S[,i]!=0)[c(which(to_mult_filt == a[1]), 
                                        which(to_mult_filt == a[2]))]
          sel_samp <- sel_samp[1]
        }
        S.copy[sel_samp, i] <- 1
      }
      return(S.copy)
    })
  }else if(poolingFilter@filter.type == "min"){
    S_all <- lapply(1:dim(X)[2], function(i){
      to_mult <- t(A.hat %*% X[,i])
      S.copy <- matrix(0, nrow = dim(S)[1], ncol = dim(S)[2])
      for(i in 1:dim(S)[2]){
        to_mult_filt <- to_mult[which(S[,i]!=0)]
        
        # Get min
        if (length(to_mult_filt) %% 2 != 0) {
          sel_samp <- which(S[,i]!=0)[which(to_mult_filt == base::min(to_mult_filt))[1]]
        } else if (length(to_mult_filt) %% 2 == 0) {
          a = sort(to_mult_filt)[c(length(to_mult_filt)/2, length(to_mult_filt)/2+1)]
          sel_samp <- which(S[,i]!=0)[c(which(to_mult_filt == a[1]), 
                                        which(to_mult_filt == a[2]))]
          sel_samp <- sel_samp[1]
        }
        S.copy[sel_samp, i] <- 1
      }
      return(S.copy)
    })
  }else if(poolingFilter@filter.type == "max"){
    S_all <- lapply(1:dim(X)[2], function(i){
      to_mult <- t(A.hat %*% X[,i])
      S.copy <- matrix(0, nrow = dim(S)[1], ncol = dim(S)[2])
      for(i in 1:dim(S)[2]){
        to_mult_filt <- to_mult[which(S[,i]!=0)]
        
        # Get max
        if (length(to_mult_filt) %% 2 != 0) {
          sel_samp <- which(S[,i]!=0)[which(to_mult_filt == base::max(to_mult_filt))[1]]
        } else if (length(to_mult_filt) %% 2 == 0) {
          a = sort(to_mult_filt)[c(length(to_mult_filt)/2, length(to_mult_filt)/2+1)]
          sel_samp <- which(S[,i]!=0)[c(which(to_mult_filt == a[1]), 
                                        which(to_mult_filt == a[2]))]
          sel_samp <- sel_samp[1]
        }
        S.copy[sel_samp, i] <- 1
      }
      return(S.copy)
    })
  }else{
    stop(paste("Pooling type", poolingFilter@filter.type, "is not defined when weights
         are applied after pooling."))
  }
  return(S_all)
}