#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param pairs A list of pairs to include in the composite model.
#' @param IntLIMresults The filtered results obtained from the IntLIM function
#' ProcessRe
#' @param modelResults A ModelResults object.
#' @param inputData An object of type IntLimData.
#' @return A final predicted value for each sample in the input data
#' @export
CompositePrediction <- function(pairs, IntLIMresults, modelResults, inputData){
  
  # Set variables for further analysis.
  covariates <- IntLIMresults@covariate.coefficients
  analyte1Vals <- inputData@analyteType1
  analyte2Vals <- inputData@analyteType2
  weights <- ComputeImportanceWeights(modelResults)
  
  # Analyte 1
  weighted_a1 <- rep(0, ncol(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_a1 <- weighted_a1 + weights[pair,] * analyte1Vals[strsplit(pair, "__")[[1]][1],]
  }
  
  # beta0
  weighted_sum_b0 <- rep(0, ncol(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b0 <- weighted_sum_b0 + weights[pair,] * rep(covariates[pair,"(Intercept)"], ncol(weights))
  }
  
  # beta1
  weighted_sum_b1 <- rep(0, ncol(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b1 <- weighted_sum_b1 + weights[pair,] * rep(covariates[pair, "a"], ncol(weights)) * 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # beta2
  weighted_sum_b2 <- rep(0, ncol(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b2 <- weighted_sum_b2 + weights[pair,] * rep(covariates[pair, "type"], ncol(weights))
  }
  
  # beta3
  weighted_sum_b3 <- rep(0, ncol(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b3 <- weighted_sum_b3 + weights[pair,] * rep(covariates[pair, "a:type"], ncol(weights))* 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # Final value.
  final_val <- (weighted_a1 - weighted_sum_b0 - weighted_sum_b1) / (weighted_sum_b2 + weighted_sum_b3)
  return(final_val)
}

#' Compute the information gain for a given prediction. The information gain
#' is computed as the difference between the entropy of the true values and the
#' conditional entropy of the true value given the prediction. Here, entropy
#' is defined using the probability that a value is greater than a selected
#' summary statistic of true values for the input data.
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param summaryStatistic The summary statistic of interest. May be "mean",
#' "median", or "midrange".
#' @export
ComputeInfoGain <- function(pred, trueVal, summaryStatistic = "mean"){
  trueVal <- trueVal[which(!is.nan(pred))]
  pred <- pred[which(!is.nan(pred))]
  originalIG <- NULL
  
  if(length(trueVal)>0 && length(pred)>0){
    middle <- mean(trueVal)
    if(summaryStatistic == "midrange"){
      middle <- min(trueVal) + (max(trueVal) - min(trueVal)) / 2
    }else if(summaryStatistic == "median"){
      middle <- median(trueVal)
    }
    
    # Compute original entropy. Note that, if the probability of one event is 0,
    # then the entropy is 0 but will be calculated as NAN.
    originalProb <- length(which(trueVal > middle)) / length(trueVal)
    originalEntropy <- -1 * (originalProb * log2(originalProb) + (1 - originalProb) * log2(1 - originalProb))
    if(is.nan(originalEntropy) && originalProb == 0){
      originalEntropy <- 0
    }
    
    # Compute entropy given prediction.
    predGreaterProb <- length(which(pred > middle)) / length(pred)
    origGreaterPredGreaterProb <- length(intersect(which(trueVal > middle), which(pred > middle))) / 
      length(pred)
    origGreaterPredGreaterTerm <- -1 * origGreaterPredGreaterProb * log2(origGreaterPredGreaterProb / predGreaterProb)
    if(is.nan(origGreaterPredGreaterTerm) == TRUE){
      origGreaterPredGreaterTerm <- 0
    }
    
    origGreaterPredLessProb <- length(intersect(which(trueVal > middle), which(pred <= middle))) / 
      length(pred)
    origGreaterPredLessTerm <- -1 * origGreaterPredLessProb * log2(origGreaterPredLessProb / (1 - predGreaterProb))
    if(is.nan(origGreaterPredLessTerm) == TRUE){
      origGreaterPredLessTerm <- 0
    }
    
    origLessPredGreaterProb <- length(intersect(which(trueVal <= middle), which(pred > middle))) /
      length(pred)
    origLessPredGreaterTerm <- -1 * origLessPredGreaterProb * log2(origLessPredGreaterProb / predGreaterProb)
    if(is.nan(origLessPredGreaterTerm) == TRUE){
      origLessPredGreaterTerm <- 0
    }
    
    origLessPredLessProb <- length(intersect(which(trueVal <= middle), which(pred <= middle))) /
      length(pred)
    origLessPredLessTerm <- -1 * origLessPredLessProb * log2(origLessPredLessProb / (1 - predGreaterProb))
    if(is.nan(origLessPredLessTerm) == TRUE){
      origLessPredLessTerm <- 0
    }
    conditionalEntropy <- origGreaterPredGreaterTerm + origGreaterPredLessTerm + origLessPredGreaterTerm + origLessPredLessTerm
    originalIG <- originalEntropy - conditionalEntropy
  }
  # If all terms are NA, then categorize this as an information loss.
  else{
    originalIG <- -1
  }
  return(originalIG)
}

#' Obtain the list of pairs in a composite model.
#' @param modelInput A ModelInput object
#' @export
ObtainCompositeSubgraphs <- function(modelInput){
  # Convert to graph.
  lineGraph <- igraph::graph_from_adjacency_matrix(modelInput@line.graph)
  
  # For each node, find its neighbors. Ensure that the neighborhood is not a 
  # subset of any of the neighbors' neighborhoods. This prevents subset graphs
  # from being returned.
  nodes <- igraph::as_ids(igraph::V(lineGraph))
  nodesToCheckQueue <- nodes
  nodesToCheckAgainst <- nodes
  
  while(length(nodesToCheckQueue) > 0){
    # Select the next node to check. Find its neighborhoods.
    node <- nodesToCheckQueue[[1]]
    neighbors <- igraph::as_ids(igraph::neighbors(lineGraph, node))
    nodeAndNeighbors <- union(neighbors, node)
    
    # For neighbors that have not yet been checked, find their neighborhoods
    # and compare to the current node's neighborhood.
    neighborsToCheck <- intersect(neighbors, nodesToCheckAgainst)
    isSubsetAny <- unlist(lapply(neighborsToCheck, function(nbr){
      isSubset <- 0
      nbrNeighbors <- igraph::as_ids(igraph::neighbors(lineGraph, nbr))
      nbrAndNeighbors <- union(nbrNeighbors, nbr)
      if(length(setdiff(nodeAndNeighbors, nbrAndNeighbors)) == 0){
        isSubset <- 1
      }
    }))
    
    # Remove the current node from the queue
    nodesToCheckQueue <- setdiff(nodesToCheckQueue, node)
    
    # If the current node is a subset of any of the other nodes, we don't
    # want to check against it anymore. Remove it from the set of nodes
    # to check against.
    if(sum(isSubsetAny) > 0){
      nodesToCheckAgainst <- setdiff(nodesToCheckAgainst, node)
    }
  }
  
  # For all remaining nodes that were not subsets, get their neighborhoods
  # and return them.
  nodeNeighborhoods <- lapply(nodesToCheckAgainst, function(node){
    neighbors <- igraph::as_ids(igraph::neighbors(lineGraph, node))
    nodeAndNeighbors <- union(neighbors, node)
  })
  return(nodeNeighborhoods)
}

#' Obtain the list of pairs in a composite model.
#' @param pairsInEachPredictor A list of pairs contained within each predictor.
#' @param importantPairs A list of all pairs that were found to be important in the
#' previous layer.
#' @export
ObtainCompositeModelsHiddenLayers <- function(pairsInEachPredictor, importantPairs){
  # Compute table of pair overlaps.
  pairCounts <- table(unlist(pairsInEachPredictor))
  predictorsContaining <- pairCounts[order(-pairCounts)]
  
  # Loop through and group composite predictors.
  belongsTo <- rep(0, length(pairsInEachPredictor))
  compositePredictors <- list()
  predictorIdx <- 1
  for(pair in names(predictorsContaining)){
    # For each pair, find the predictors containing it.
    listContainsPair <- unlist(lapply(1:length(pairsInEachPredictor), function(i){
      composite <- pairsInEachPredictor[[i]]
      containsPair <- FALSE
      if(pair %in% composite){
        containsPair <- TRUE
      }
      return(containsPair)
    }))
    # If the predictors containing it have not been assigned, assign them.
    if(max(belongsTo[which(listContainsPair == TRUE)]) == 0){
      belongsTo[which(listContainsPair == TRUE)] <- predictorIdx
      predictorIdx <- predictorIdx + 1
    }
    # If some of the predictors have already been assigned, merge with them.
    else{
      belongsTo[which(listContainsPair == TRUE)] <- max(belongsTo[which(listContainsPair == TRUE)])
    }
  }
  pairList <- lapply(1:length(unique(belongsTo)), function(i){
    return(intersect(names(predictorsContaining)[which(belongsTo == i)], importantPairs))
  })
  return(pairList)
}

#' Given all predictors, find only the important predictors.
#' @param pairs A list of pairs to include in the composite model.
#' @param IntLIMresults The filtered results obtained from the IntLIM function
#' ProcessRe
#' @param modelResults A ModelResults object.
#' @param inputData An object of type IntLimData.
#' @param verbose Whether or not to print out each step.
#' @return A list of informative pairs
#' @export
PrunePredictors <- function(compositeSubgraphs, IntLIMresults, modelResults, inputData, verbose = FALSE){
  prunedSubgraphs <- lapply(1:length(compositeSubgraphs), function(i){
    # Print initial graph.
    if(verbose == TRUE){
      print(unlist(compositeSubgraphs[[i]]))
    }
    
    # Initialize the predictor to include all pairs in the composite subgraph.
    importantPairs <- compositeSubgraphs[[i]]
    compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = importantPairs, 
                                                                     IntLIMresults = IntLIMresults,
                                                                     modelResults = modelResults,
                                                                     inputData = inputData)
    infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel), 
                                                           trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
    removedLastTime <- FALSE
    compositeModelFull <- compositeModel
    infoGainFull <- infoGain
    
    # Sequentially test removal of each pair.
    for(pair in compositeSubgraphs[[i]]){
      # Save the information gain for the full model.
      if(removedLastTime == TRUE){
        compositeModelFull <- compositeModel
        infoGainFull <- infoGain
      }
      
      # Compute the information gain for the new model.
      pairsToInclude <- setdiff(importantPairs, pair)
      if(length(pairsToInclude) > 0){
        compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = pairsToInclude,
                                                                         IntLIMresults = IntLIMresults,
                                                                         modelResults = modelResults,
                                                                         inputData = inputData)
        infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel),
                                                               trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
      }
      
      # If the information gain for the new model is greater than for the full model, remove the pair.
      removedLastTime <- FALSE
      if(infoGain >= infoGainFull && length(pairsToInclude) > 0){
        # Print the increase in information gain.
        if(verbose == TRUE){
          print(paste(list("Removed pair. Information gain after removing", pair, "is", infoGain), collapse = " "))
        }
        
        importantPairs <- pairsToInclude
        removedLastTime <- TRUE
        
        tryCatch({
          # Plot the new line graph.
          newPredictions <- t(matrix(rep(compositeModel, nrow(modelResults@model.input@node.wise.prediction)),
                                     ncol = nrow(modelResults@model.input@node.wise.prediction)))
          rownames(newPredictions) <- rownames(modelResults@model.input@node.wise.prediction)
          colnames(newPredictions) <- colnames(modelResults@model.input@node.wise.prediction)
          modelResultsNew <- modelResults
          modelResultsNew@model.input@node.wise.prediction <- newPredictions
          if(verbose == TRUE){
            MultiOmicsGraphPrediction::PlotLineGraph(modelResults = modelResultsNew, subset = pairsToInclude,
                                                     sampSubset = "UACC.62", weights = TRUE, stype = "drug5FU")
          }
        }, warning=function(cond){
          if(verbose == TRUE){
            plot.new()
          }
        })
        
      }else{
        # Print the decrease in information gain.
        if(verbose == TRUE){
          print(paste(list("Kept pair. Information gain after removing", pair, "is", infoGain), collapse = " "))
        }
      }
    }
    return(importantPairs)
  })
  return(prunedSubgraphs)
}

#' #' Given multiple composite predictors, prune the predictors that are not needed.
#' #' @param pairs A list of pairs to include in the composite model.
#' #' @param IntLIMresults The filtered results obtained from the IntLIM function
#' #' @param modelResults A ModelResults object.
#' #' @param inputData An object of type IntLimData.
#' #' @return A list of informative pairs
#' #' @export
#' PrunePredictorsHiddenLayers <- function(pairs, IntLIMresults, modelResults, inputData){
#'   importantPredictors <- lapply(pairs, function(pairsPred){
#'     pred1 <- CompositePrediction(pairs = pairs, 
#'                                  IntLIMresults = IntLIMresults, 
#'                                  modelResults = modelResults,
#'                                  inputData = inputData)
#'     
#'     # Prune the pairs.
#'     removePairPred <- CompositePredictionLOO(pairs = pairs, 
#'                                              IntLIMresults = IntLIMresults, 
#'                                              modelResults = modelResults,
#'                                              inputData = inputData)
#'     infoGainOriginal <- ComputeInfoGain(pred = unlist(pred1), 
#'                                         trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
#'     infoGainsLOO <- lapply(removePairPred, function(p){
#'       return(ComputeInfoGain(pred = unlist(p), 
#'                              trueVal = inputData@sampleMetaData[,IntLIMresults@stype]))
#'     })
#'     informativePairs <- unlist(pairsPred[which(infoGainsLOO < infoGainOriginal)])
#'     if(length(which(infoGainsLOO < infoGainOriginal))==0){
#'       informativePairs <- unlist(pairsPred[which(infoGainsLOO == infoGainOriginal)])[1]
#'     }
#'     
#'     # For each pair that was pruned, see if it functions better than the full
#'     # predictor as an "orphan" pair.
#'     orphans <- setdiff(unlist(pairs), informativePairs)
#'     importantOrphans <- unlist(lapply(1:length(orphans), function(orphan){
#'       orphanPred <- pred[,orphan]
#'       orphanInfoGain <- ComputeInfoGain(pred = orphanPred, 
#'                                         trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
#'       retval <- NULL
#'       if(orphanInfoGain > infoGainOriginal){
#'         retval <- orphans[orphan]
#'       }
#'       return(retval)
#'     }))
#'     return(list(modelPairs = informativePairs, orphans = importantOrphans))
#'   })
#' }

#' Compute the weight of each predictor given the weights of different
#' importance metrics.
#' @param modelResults A ModelResults object.
#' @return A list of informative pairs
#' @export
ComputeImportanceWeights <- function(modelResults){
  importanceWeights <- modelResults@current.importance.weights
  importanceWeightMatrix <- t(matrix(rep(importanceWeights, nrow(modelResults@importance[[1]])), 
                                     ncol = nrow(modelResults@importance[[1]]),
                                     nrow = ncol(modelResults@importance[[1]])))
  weights <- lapply(1:length(modelResults@importance), function(i){
    imp <- modelResults@importance[[i]] * importanceWeightMatrix
    nm <- names(modelResults@importance)[i]
    sum_over <- as.data.frame(rowSums(imp))
    rownames(sum_over) <- rownames(imp)
    colnames(sum_over) <- nm
    return(sum_over)
  })
  weights_all <- do.call(cbind,weights)
  wt <- as.matrix(weights_all)
  return(wt)
}
