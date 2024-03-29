#' Find edges that share nodes and add them to a data frame.
#' @param modelInputs An object of the ModelInput class.
#' @param iterations Maximum number of iterations. Default is 1,000.
#' @param convergenceCutoff Cutoff for convergence. Default is 0.001.
#' @param learningRate Learning rate to use during training. Default is 0.2
#' @param activationType Activation function. May be "softmax", "sigmoid", 
#' "tanh", or "none". Default is "none", meaning stype.class is numeric.
#' @param optimizationType Type of optimization. May be "SGD", "momentum",
#' "adagrad", or "adam". Default is "SGD".
#' @param initialMetaFeatureWeights Initial weights for model meta-features. Default
#' is 0, which results in each meta-feature being given equal weight.
#' @return An object of the ModelResults class. This object will later be filled in
#' with the weights and errors from each iteration as well as the initial
#' settings and data and the final result.
#' @export
InitializeGraphLearningModel <- function(modelInputs,
                                         iterations = 1000,
                                         convergenceCutoff = 0.001,
                                         learningRate = 0.2,
                                         activationType = "none", 
                                         optimizationType = "SGD",
                                         initialMetaFeatureWeights = 0){
  # Initialize metafeature weights.
  weights_count <- length(modelInputs@metaFeatures)
  wt_name <- names(modelInputs@metaFeatures)
  
  # Initialize data frame with maximum number of iterations.
  tracking.frame <- as.data.frame(matrix(-1, nrow = iterations, 
                                         ncol = 2 + (2 * weights_count)))
  tracking.frame.cnames <- c("Iteration", "Error")
  tracking.frame.cnames <- c(tracking.frame.cnames, paste("Weight", wt_name, sep = "_"))
  tracking.frame.cnames <- c(tracking.frame.cnames, paste("Gradient", wt_name, sep = "_"))
  colnames(tracking.frame) <- tracking.frame.cnames
  tracking.frame$Error[1] <- .Machine$double.xmax
  tracking.frame$Iteration[1] <- 0

  # Initialize weights with uniform distribution.
  max_phen <- max(modelInputs@true.phenotypes)
  num_nodes <- nrow(modelInputs@edge.wise.prediction)
  weights <- as.matrix(rep(1 / weights_count, weights_count))
  if(initialMetaFeatureWeights != 0){
    weights <- initialMetaFeatureWeights
  }
  names(weights) <- names(modelInputs@metaFeatures)
  
  tracking.frame[1,3:(2+weights_count)] <- weights

  # Initialize and return results.
  newModelResults <- methods::new("ModelResults", model.input=modelInputs,
                        iteration.tracking=tracking.frame, max.iterations=iterations,
                        convergence.cutoff=convergenceCutoff, learning.rate=learningRate,
                        previous.metaFeature.weights=as.matrix(rep(0,length(weights))), 
                        current.metaFeature.weights=as.matrix(weights),
                        current.gradient=as.matrix(rep(-1,length(weights))),
                        previous.momentum=as.matrix(rep(0,length(weights))),
                        previous.update.vector=as.matrix(rep(0,length(weights))),
                        sum.square.gradients=as.matrix(rep(0,length(weights))),
                        current.iteration=0, activation.type=activationType,
                        optimization.type=optimizationType,
                        pairs = "")
  return(newModelResults)
}

#' 
#' Format the input for graph-based learning. This input consists of:
#' - The Laplacian of a line graph built from the co-regulation graphs, where 
#' each node corresponds to a pair of analytes.
#' - A prediction value for each node of the line graph, for each sample X.
#' - The true prediction values Y for each sample X.
#' - The co-regulation graph.
#' - The line graph built from the co-regulation graphs.
#' - The input data (analyte levels, covariates, and phenotype).
#' - The model results from IntLIM for each predictor.
#' - The metafeatures for each model from IntLIM.
#' - The phenotype / outcome column name in the input data.
#' - The type of the outcome / phenotype ("numeric" or "categorical").
#' - The analyte type used as the outcome variable for each model (1 or 2).
#' - The analyte type used as the independent variable for each model (1 or 2).
#' @slot metaFeatures A list of calculated meta-features for each sample
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
#' @param edgeTypeList List containing one or more of the following to include
#' in the line graph:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
#' @param verbose Whether to print the number of predictions replaced in each sample.
#' TRUE or FALSE. Default is FALSE.
#' @param outcome The outcome used in the IntLIM models.
#' @param independent.var.type The independent variable type used in the IntLIM models.
#' @return An object of the ModelInput class. This object will contain the input data,
#' graph information, and individual predictor properties and their metafeatures.
#' @export
FormatInput <- function(predictionGraphs, coregulationGraph, metaFeatures, modelProperties,
                        inputData, stype.class, edgeTypeList, stype, verbose = FALSE, covariates = list(),
                        outcome = 2, independent.var.type = 2){
  
  # Extract edge-wise predictions.
  predictions_by_edge <- lapply(names(predictionGraphs), function(sampName){
    df_predictions <- igraph::as_data_frame(predictionGraphs[[sampName]])
    edge_names <- paste(make.names(df_predictions$from), make.names(df_predictions$to),
                        sep = "__")
    df_predictions_new <- data.frame(Edge = edge_names, Weight = df_predictions$weight)
    return(df_predictions_new)
  })
  names(predictions_by_edge) <- names(predictionGraphs)
  predicted_weights_only <- lapply(predictions_by_edge, function(pred){
    return(pred$Weight)
  })
  predictions_flattened <- t(data.frame(predicted_weights_only))
  colnames(predictions_flattened) <- predictions_by_edge[[1]]$Edge
  
  # Obtain the predictions.
  Y <- inputData@sampleMetaData[,stype]
  if(stype.class == "factor"){
    Y <- as.numeric(Y)-1
  }
  names(Y) <- names(predictions_by_edge)
  
  # Create a ModelInput object and return it.
  newModelInput <- methods::new("ModelInput", edge.wise.prediction=predictions_flattened,
                                true.phenotypes=Y, coregulation.graph=igraph::get.adjacency(coregulationGraph, sparse = FALSE), 
                                input.data = inputData, model.properties = modelProperties,
                                metaFeatures = metaFeatures, stype = stype, covariates = covariates, stype.class = stype.class,
                                outcome = outcome, independent.var.type = independent.var.type)
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
#' @return An adjacency matrix for the line graph.
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
#' @return A data frame of all nodes (analyte pairs / predictors) that
#' share an analyte, with the columns "to" and "from".
FindEdgesSharingNodes <- function(predictionsByEdge, graphWithPredictions, nodeType1,
                                  nodeType2){
  # Convert predictions to data frame.
  graph_df <- igraph::as_data_frame(graphWithPredictions)

  # Find shared nodes.
  nodes <- unique(graph_df[,nodeType1])
  to_shared <- lapply(1:length(nodes), function(i){
    
    # Find all line graph nodes starting with or ending with the analyte in question.
    node <- nodes[i]
    combs <- 0
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
  str(line_graph_df)
  
  return(line_graph_df)
}

#' Optimize the combination of predictors by metafeatures alone (in other words,
#' exclude pooling and combine in a single layer using a linear combination).
#' @param modelResults An object of the ModelResults class.
#' @param verbose Whether to print results as you run the model.
#' @param modelRetention Strategy for model retention. "stringent" (the default)
#' retains only models that improve the prediction score. "lenient" also retains models that
#' neither improve nor reduce the prediction score.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @return A ModelResults object, with all of the tracking information from each
#' iteration filled in.
#' @export
OptimizeMetaFeatureCombo <- function(modelResults, verbose = TRUE, modelRetention = "stringent",
                                    useCutoff = FALSE){
  
  # Calculate prediction cutoffs.
  predictionCutoffs <- CalculatePredictionCutoffs(modelResults = modelResults)
  
  # Start the first iteration and calculate a dummy weight delta.
  modelResults@current.iteration <- 1
  weight.delta <- sqrt(sum((modelResults@current.metaFeature.weights - 
                              modelResults@previous.metaFeature.weights)^2))
  
  # Get the initial pruned models.
  pruningMethod <- "error.t.test"
  pairsPredAll <- MultiOmicsGraphPrediction::ObtainSubgraphNeighborhoods(modelInput = modelResults@model.input, percentOverlapCutoff = 50)
  prunedModels <- MultiOmicsGraphPrediction::DoSignificancePropagation(pairs = pairsPredAll, modelResults = modelResults,
                                                                          verbose = FALSE,
                                                                       pruningMethod = pruningMethod,
                                                                       modelRetention = modelRetention,
                                                                       useCutoff = useCutoff,
                                                                       minCutoff = predictionCutoffs$min,
                                                                       maxCutoff = predictionCutoffs$max)
  
  # Set initial error.
  Y.pred <- DoPrediction(modelResults = modelResults, prunedModels = prunedModels,
                         useCutoff = useCutoff,
                         minCutoff = predictionCutoffs$min,
                         maxCutoff = predictionCutoffs$max)
  if(modelResults@model.input@stype.class == "categorical"){
    modelResults@iteration.tracking$Error[1] <- 
      ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
  }else{
    modelResults@iteration.tracking$Error[1] <- 
      ComputeRMSE(modelResults@model.input@true.phenotypes, Y.pred)
  }
  print(paste("Initial error is", modelResults@iteration.tracking$Error[1]))
  
  # Repeat the training process for all iterations, until the maximum is reached
  # or until convergence.
  sequential_convergence_count <- 0
  sequential_count_limit <- 5
  while(modelResults@current.iteration <= modelResults@max.iterations
        && (sequential_convergence_count < sequential_count_limit)){
    # For stochastic training, permute the samples, then compute the gradient one
    # sample at a time.
    perm_samples <- sample(1:nrow(modelResults@model.input@node.wise.prediction),
                          nrow(modelResults@model.input@node.wise.prediction))

    # Initialize previous weight vector.
    modelResults@previous.metaFeature.weights <- modelResults@current.metaFeature.weights
    
    for(i in perm_samples){
      # Extract data for each sample.
      newModelResults <- modelResults
      newModelResults@model.input@node.wise.prediction <- 
        as.matrix(modelResults@model.input@node.wise.prediction[i,])
      newModelResults@model.input@true.phenotypes <- 
        modelResults@model.input@true.phenotypes[i]
      newModelResults@model.input@input.data@analyteType1 <- as.matrix(modelResults@model.input@input.data@analyteType1[,i])
      newModelResults@model.input@input.data@analyteType2 <- as.matrix(modelResults@model.input@input.data@analyteType2[,i])
      newModelResults@model.input@input.data@sampleMetaData <- as.data.frame(modelResults@model.input@input.data@sampleMetaData[i,])
      metaFeaturesSamp <- lapply(1:length(modelResults@model.input@metaFeatures), function(imp){
        df <- t(as.data.frame(modelResults@model.input@metaFeatures[[imp]][i,]))
        rownames(df) <- rownames(modelResults@model.input@node.wise.prediction)[i]
        colnames(df) <- colnames(modelResults@model.input@node.wise.prediction)
        return(df)
      })
      names(metaFeaturesSamp) <- names(modelResults@model.input@metaFeatures)
      newModelResults@model.input@metaFeatures <- metaFeaturesSamp
      
      # Do training iteration.
      newModelResults <- DoSingleTrainingIteration(modelResults = newModelResults,
                                                   prunedModels = prunedModels,
                                                   useCutoff = useCutoff,
                                                   minCutoff = predictionCutoffs$min,
                                                   maxCutoff = predictionCutoffs$max)
      
      # Update weights and gradient in the model results according to the
      # results of this sample.
      modelResults@current.metaFeature.weights <- newModelResults@current.metaFeature.weights
      modelResults@previous.metaFeature.weights <- newModelResults@previous.metaFeature.weights
      modelResults@current.gradient <- newModelResults@current.gradient
      modelResults@previous.momentum <- newModelResults@previous.momentum
      modelResults@previous.update.vector <- newModelResults@previous.update.vector
      modelResults@iteration.tracking <- newModelResults@iteration.tracking
    }
    # Fill in the prediction for error calculation.
    Y.pred <- DoPrediction(modelResults = modelResults, prunedModels = prunedModels,
                           useCutoff = useCutoff,
                           minCutoff = predictionCutoffs$min,
                           maxCutoff = predictionCutoffs$max)

    # Compute the prediction error over all samples.
    modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1] <- modelResults@current.iteration
    if(modelResults@model.input@stype.class == "categorical"){
      modelResults@iteration.tracking$Error[modelResults@current.iteration+1] <- 
        ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
    }else{
      modelResults@iteration.tracking$Error[modelResults@current.iteration+1] <- 
        ComputeRMSE(modelResults@model.input@true.phenotypes, Y.pred)
    }
    currentError <- modelResults@iteration.tracking$Error[modelResults@current.iteration+1]
    
    # Get the new pruned models.
    prunedModels <- MultiOmicsGraphPrediction::DoSignificancePropagation(pairs = pairsPredAll, modelResults = modelResults,
                                                                        verbose = FALSE,
                                                                        modelRetention = modelRetention,
                                                                        useCutoff = useCutoff,
                                                                        minCutoff = predictionCutoffs$min,
                                                                        maxCutoff = predictionCutoffs$max)
    modelResults@pairs <- unlist(prunedModels)
    prunedModelSizes <- lapply(prunedModels, function(model){return(length(model))})
    
    # Print the weight delta and error.
    weight.delta <- sqrt(sum((modelResults@current.metaFeature.weights - modelResults@previous.metaFeature.weights)^2))
    if(modelResults@current.iteration %% 1 == 0 && verbose == TRUE){
      print(paste("iteration", modelResults@current.iteration, ": weight delta is", weight.delta,
                  "and error is", paste0(currentError, ". Final subgraph has ", prunedModelSizes, " nodes.")))
      sortedNodes <- sort(unlist(prunedModels))
      if(prunedModelSizes > 5){
        sortedNodes <- sortedNodes[1:5]
      }
      print(paste0("Nodes include: ", paste(sortedNodes, collapse = ", ")))
    }
    
    # Update the iteration.
    modelResults@current.iteration <- modelResults@current.iteration + 1
    modelResults@pairs <- unlist(prunedModels)
    
    # Increment the number of convergent iterations if applicable.
    if(weight.delta < modelResults@convergence.cutoff){
      sequential_convergence_count = sequential_convergence_count + 1
    }else{
      sequential_convergence_count = 0
    }
  }
  # If we exited before the maximum number of iterations, remove the rest of the
  # tracking data.
  if((modelResults@current.iteration-1) < modelResults@max.iterations){
    modelResults@iteration.tracking <- modelResults@iteration.tracking[1:(modelResults@current.iteration-1),]
  }
  return(modelResults)
}

#' Compute classification error.
#' @param true.Y The true phenotype of each sample.
#' @param pred.Y The predicted phenotype of each sample.
#' @return A vector of classification errors.
#' @export
ComputeClassificationError <- function(true.Y, pred.Y){
  # Find false and true positives and negatives.
  FP <- length(intersect(which(true.Y == 1), which(pred.Y == 2)))
  TP <- length(intersect(which(true.Y == 2), which(pred.Y == 2)))
  FN <- length(intersect(which(true.Y == 2), which(pred.Y == 1)))
  TN <- length(intersect(which(true.Y == 1), which(pred.Y == 1)))

  # Compute error and return.
  error <- (FP + FN) / (FP + FN + TP + TN)
  return(error)
}

#' Compute the normalized root mean squared error.
#' @param true.Y The true phenotype of each sample.
#' @param pred.Y The predicted phenotype of each sample.
#' @return A vector of RMSE.
#' @export
ComputeRMSE <- function(true.Y, pred.Y){
  RMSD <- sqrt(sum((true.Y - pred.Y)^2) / length(true.Y))
  return(RMSD)
}