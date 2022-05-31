#' 
#' Format the input for graph-based learning. This input consists of:
#' 1. The Laplacian of a line graph built from the co-regulation graphs, where 
#' each node corresponds to a pair of analytes.
#' 2. A prediction value for each node of the line graph, for each sample X.
#' 3. The true prediction values Y for each sample X.
#' @slot importance A list of calculated importance metric data frames for each sample
#' @slot modelProperties A data frame that includes model information, i.e. R^2,
#' interaction term p-value, and coefficients.
#' @slot inputData An IntLIMData object that includes slots for the sample data,
#' the analyte data, and the analyte meta data.
#' @param predictionGraphs A list of igraph objects, each of which includes
#' predictions for each edge.
#' @param coregulationGraph An igraph object containing the coregulation graph.
#' @param stype.class The class of the outcome ("numeric" or "categorical")
#' @param stype The phenotype or outcome of interest
#' @param covariates A list of covariates to include in the model. These will be in the
#' sampleMetaData slot of the inputData variable.
#' @param predictorBounds Upper and lower bounds for predictors. Predictors outside
#' of these bounds for a sample will not be included in the final prediction for that sample. 
#' Defaults to c(NA,NA), meaning that no predictors will be removed.
#' @param edgeTypeList List containing one or more of the following to include
#' in the line graph:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
#' @param verbose Whether to print the number of predictions replaced in each sample.
#' TRUE or FALSE. Default is FALSE.
#' @export
FormatInput <- function(predictionGraphs, coregulationGraph, importance, modelProperties,
                        inputData, stype.class, edgeTypeList, stype, verbose = FALSE, covariates = list(),
                        predictorBounds = c(NA,NA)){
  
  # Extract edge-wise predictions.
  predictions_by_node <- lapply(names(predictionGraphs), function(sampName){
    df_predictions <- igraph::as_data_frame(predictionGraphs[[sampName]])
    node_names <- paste(make.names(df_predictions$from), make.names(df_predictions$to),
                        sep = "__")
    df_predictions_new <- data.frame(Node = node_names, Weight = df_predictions$weight)
    return(df_predictions_new)
  })
  names(predictions_by_node) <- names(predictionGraphs)
  predicted_weights_only <- lapply(predictions_by_node, function(pred){
    return(pred$Weight)
  })
  predictions_flattened <- t(data.frame(predicted_weights_only))
  colnames(predictions_flattened) <- predictions_by_node[[1]]$Node
  
  # Convert co-regulation graph into a line graph. Return the adjacency matrix.
  # If edges were not connected by nodes in the original graph, they may be
  # removed from the line graph. Remove these from the predictions_by_node df.
  A <- CreateLineGraph(predictionsByEdge = predictions_by_node[[1]],
                       graphWithPredictions = predictionGraphs[[1]],
                       edgeTypeList = edgeTypeList)
  predictions_flattened <- predictions_flattened[,colnames(A)]
  predictions_flattened_orig <- predictions_flattened
  
  # Add self-loops.
  A_tilde <- as.matrix(A)
  diag(A_tilde) <- 1
  
  # Extract the diagonal (in degree) and raise to the negative half power.
  diags1 <- colSums(A_tilde)
  diags_neg_half <- 1 / sqrt(diags1)
  D_tilde_neg_half1 <- t(matrix(rep(diags_neg_half,length(diags_neg_half)),
                                nrow = length(diags_neg_half)))
  A_hat1 <- D_tilde_neg_half1 * A_tilde
  rm(D_tilde_neg_half1)
  
  # Extract the diagonal (out degree) and raise to the negative half power.
  diags2 <- rowSums(A_tilde)
  rm(A_tilde)
  diags_neg_half <- 1 / sqrt(diags2)
  D_tilde_neg_half2 <- t(matrix(rep(diags_neg_half,length(diags_neg_half)),
                                nrow = length(diags_neg_half)))
  
  # Obtain the final matrix. Note that we modify the matrix multiplication
  # problem to obtain an elementwise multiplication problem
  # because it speeds up computation.
  A_hat <- A_hat1 * D_tilde_neg_half2
  
  # Obtain the predictions.
  Y <- inputData@sampleMetaData[,stype]
  if(stype.class == "factor"){
    Y <- as.numeric(Y)-1
  }
  names(Y) <- names(predictions_by_node)
  
  # Create a binary predictor mask. The mask should include a "1" for every predictor
  # within the specified bounds for a given sample and a "0" for every predictor out of bounds.
  predictor_mask <- matrix(rep(1, times = ncol(predictions_flattened) * nrow(predictions_flattened)),
                           ncol = ncol(predictions_flattened))
  rownames(predictor_mask) <- rownames(predictions_flattened)
  colnames(predictor_mask) <- colnames(predictions_flattened)
  if(!is.na(predictorBounds[1]) && !is.na(predictorBounds[2])){
    predictor_mask[which(predictions_flattened < predictorBounds[1])] <- 0
    predictor_mask[which(predictions_flattened > predictorBounds[2])] <- 0
  }
  
  # Create a ModelInput object and return it.
  newModelInput <- methods::new("ModelInput", A.hat=A_hat, node.wise.prediction=predictions_flattened,
                                true.phenotypes=Y, outcome.type=stype.class, 
                                coregulation.graph=igraph::get.adjacency(coregulationGraph, sparse = FALSE), 
                                line.graph=as.matrix(A), input.data = inputData, model.properties = modelProperties,
                                importance = importance, stype = stype, covariates = covariates, stype.class = stype.class,
                                mask = predictor_mask)
  return(newModelInput)
}

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
    combs <- expand.grid(set_with_2, set_with_1)
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
#' @param prunedModels The models that remain after pruning. There should only
#' be one model at the end.
DoPrediction <- function(modelResults, prunedModels){
  
  # Obtain the final prediction.
  predictions <- CompositePrediction(pairs = prunedModels, modelResults = modelResults)
  Y.pred <- as.data.frame(predictions)
  colnames(Y.pred) <- 1
  rownames(Y.pred) <- names(predictions)

  # Use activation function if output is of a character type. Note that all character
  # types are converted into factors, and since only binary factors are accepted by
  # the package, the values will be 1 (for the alphanumerically lowest level) and 2
  # (for the alphanumerically highest level).
  if(modelResults@model.input@outcome.type == "factor"){
    if(modelResults@activation.type == "softmax"){
      Y.pred <- round(SoftmaxWithCorrection(Y.pred))
    }else if(modelResults@activation.type == "tanh"){
      Y.pred <- round(TanhWithCorrection(Y.pred))
    }else{
      Y.pred <- round(SigmoidWithCorrection(Y.pred))
    }
  }
  
  # Return the prediction.
  return(Y.pred)
}

#' Train the graph learning model, using the specifications in the ModelResults
#' class and storing the results in the ModelResults class.
#' @param modelResults An object of the ModelResults class.
#' @param prunedModels The models that remain after pruning.
DoSingleTrainingIteration <- function(modelResults, prunedModels){
  # Predict Y.
  Y.pred <- DoPrediction(modelResults, prunedModels)
  modelResults@outcome.prediction <- as.numeric(Y.pred)

  # Backpropagate and calculate the error.
  error <- modelResults@iteration.tracking$Error[modelResults@current.iteration-1]
  modelResults <- Backpropagate(modelResults = modelResults,
                                prunedModels = prunedModels,
                                Y.pred = Y.pred)
  if(modelResults@model.input@outcome.type == "categorical"){
    modelResults@iteration.tracking$Error[iteration+1] <- 
      ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
  }else{
    modelResults@iteration.tracking$Error[modelResults@current.iteration+1] <- 
      ComputeRMSE(modelResults@model.input@true.phenotypes, Y.pred)
  }
  
  # Modify the model results and return.
  return(modelResults)
}