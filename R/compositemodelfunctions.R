#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param pairs A list of pairs to include in the composite model.
#' @param IntLIMresults The filtered results obtained from the IntLIM function
#' ProcessRe
#' @param modelResults A ModelResults object.
#' @param inputData An object of type IntLimData.
#' @param covar Covariates to be used in the model.
#' @return A final predicted value for each sample in the input data
#' @export
CompositePrediction <- function(pairs, IntLIMresults, modelResults, inputData, covar = NULL){
  
  # Set variables for further analysis.
  covariates <- IntLIMresults@covariate.coefficients
  analyte1Vals <- inputData@analyteType1
  analyte2Vals <- inputData@analyteType2
  covariateVals <- inputData@sampleMetaData
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
  
  # covariates
  weighted_sum_covars <- rep(0, ncol(weights))
  if(!is.null(covar)){
    weighted_sum_each <- lapply(covar, function(c){
      weighted_sum <- rep(0, ncol(weights))
      for(i in 1:length(pairs)){
        pair <- pairs[[i]]
        weighted_sum <- weighted_sum + weights[pair,] * rep(covariates[pair, c], ncol(weights)) *
          covariateVals[,c]
      }
      return(weighted_sum)
    })
    weighted_sum_covars <- Reduce('+', weighted_sum_each)
  }
  
  # Final value.
  final_val <- (weighted_a1 - weighted_sum_b0 - weighted_sum_b1 - weighted_sum_covars) / (weighted_sum_b2 + weighted_sum_b3)
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
ObtainCompositeModelsFirstLayer <- function(pairsInEachPredictor, importantPairs){
  # Compute table of pair overlaps.
  pairCounts <- table(unlist(pairsInEachPredictor))
  predictorsContaining <- pairCounts[order(-pairCounts)]
  
  # Loop through and group composite predictors.
  belongsTo <- rep(0, length(pairsInEachPredictor))
  compositePredictors <- importantPairs
  compositeFull <- pairsInEachPredictor
  pairMapping <- data.frame(from = 1:length(importantPairs), to = 1:length(importantPairs))
  for(pair in names(predictorsContaining[which(predictorsContaining > 1)])){
    # For each pair, find the predictors containing it.
    listContainsPair <- unlist(lapply(1:length(compositePredictors), function(i){
      composite <- compositePredictors[[i]]
      containsPair <- FALSE
      if(pair %in% composite){
        containsPair <- TRUE
      }
      return(containsPair)
    }))
    whichContainsPair <- which(listContainsPair == TRUE)
    if(length(whichContainsPair) > 0){
      
      # Create a new combined predictor.
      newPredictor <- unique(unlist(compositePredictors[whichContainsPair]))
      newFull <- unique(unlist(compositeFull[whichContainsPair]))
      
      # Replace old predictor information with new.
      firstPiece <- min(whichContainsPair)
      compositePredictors[[firstPiece]] <- newPredictor
      compositeFull[[firstPiece]] <- newFull
      toKeep <- sort(union(firstPiece, setdiff(1:length(compositePredictors), whichContainsPair)))
      compositePredictors <- compositePredictors[toKeep]
      compositeFull <- compositeFull[toKeep]
      
      # Update the mappings.
      pairMapping[which(pairMapping$to %in% whichContainsPair), "to"] <- min(whichContainsPair)
      indicesToUpdate <- which(pairMapping$to > max(whichContainsPair))
      pairMapping[indicesToUpdate, "to"] <- pairMapping[indicesToUpdate, "to"] - (length(whichContainsPair) - 1)
    }
  }
  return(list(compositeModels = compositePredictors, expandedCompositeModels = compositeFull, mapping = pairMapping))
}

#' Given all predictors, find only the important predictors.
#' @param compositeSubgraphs A list of pairs to include in the composite model.
#' @param IntLIMresults The filtered results obtained from the IntLIM function
#' ProcessRe
#' @param modelResults A ModelResults object.
#' @param inputData An object of type IntLimData.
#' @param verbose Whether or not to print out each step.
#' @param covar Covariates to be used in the model.
#' @return A list of informative pairs
#' @export
PrunePredictors <- function(compositeSubgraphs, IntLIMresults, modelResults, inputData, covar = NULL, verbose = FALSE){
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
                                                                     covar = covar,
                                                                     inputData = inputData)
    infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel), 
                                                           trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
    print(paste(list("Original information gain is", infoGain), collapse = " "))
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
                                                                         covar = covar,
                                                                         inputData = inputData)
        infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel),
                                                               trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
      }
      
      # If the information gain for the new model is greater than for the full model, remove the pair.
      removedLastTime <- FALSE
      if(infoGain >= infoGainFull && length(pairsToInclude) > 0){
        # Print the increase in information gain.
        if(verbose == TRUE){
          print(paste(list("Removed pair. Information gain after removing", pair, "is", infoGain - infoGainFull), collapse = " "))
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
          print(paste(list("Kept pair. Information loss after removing", pair, "is", infoGainFull - infoGain), collapse = " "))
        }
      }
    }
    return(importantPairs)
  })
  return(prunedSubgraphs)
}

#' Given multiple composite predictors, prune the predictors that are not needed.
#' @param compositeSubgraphs A list of pairs to include in the composite model.
#' @param previousModels A list of the previous models that were consolidated.
#' @param IntLIMresults The filtered results obtained from the IntLIM function
#' ProcessRe
#' @param modelResults A ModelResults object.
#' @param inputData An object of type IntLimData.
#' @param covar Covariates to be used in the model.
#' @param verbose Whether or not to print out each step.
#' @return A list of informative pairs
#' @export
PrunePredictorsHiddenLayers <- function(compositeSubgraphs, previousModels, IntLIMresults, modelResults, inputData, covar = NULL, verbose = FALSE){
  # Extract relevant information.
  pairs <- compositeSubgraphs$compositeModels
  mapping <- compositeSubgraphs$mapping
  
  prunedSubgraphs <- lapply(1:length(pairs), function(i){
    # Print initial graph.
    if(verbose == TRUE){
      print(sort(unlist(pairs[[i]])))
    }
    
    # Initialize the predictor to include all pairs in the composite subgraph.
    importantPairs <- pairs[[i]]
    compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = importantPairs, 
                                                                     IntLIMresults = IntLIMresults,
                                                                     modelResults = modelResults,
                                                                     covar = covar,
                                                                     inputData = inputData)
    infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel), 
                                                           trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
    print(paste(list("Original information gain is", infoGain), collapse = " "))
    removedLastTime <- FALSE
    compositeModelFull <- compositeModel
    infoGainFull <- infoGain
    
    # Figure out which of the previous models mapped to this one.
    previousModelsMapped <- mapping$from[which(mapping$to == i)]
    
    # Sequentially test removal of each model.
    for(model in previousModelsMapped){
      
      # Extract the list of pairs to keep.
      modelsToUse <- setdiff(previousModelsMapped, model)
      pairsToInclude <- Reduce(union, previousModels[modelsToUse])

      # Save the information gain for the full model.
      if(removedLastTime == TRUE){
        compositeModelFull <- compositeModel
        infoGainFull <- infoGain
      }
      
      # Compute the information gain for the new model.
      if(length(pairsToInclude) > 0){
        compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = pairsToInclude,
                                                                         IntLIMresults = IntLIMresults,
                                                                         modelResults = modelResults,
                                                                         covar = covar,
                                                                         inputData = inputData)
        infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel),
                                                               trueVal = inputData@sampleMetaData[,IntLIMresults@stype])
      }
      
      # If the information gain for the new model is greater than for the full model, remove the pair.
      removedLastTime <- FALSE
      if(infoGain >= infoGainFull && length(pairsToInclude) > 0){
        # Print the increase in information gain.
        if(verbose == TRUE){
          print(paste(list("Removed model. Information gain after removing", model, "is", infoGain - infoGainFull), collapse = " "))
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
          print(paste(list("Kept model. Information loss after removing", model, "is", infoGainFull - infoGain), collapse = " "))
        }
      }
    }
    return(importantPairs)
  })
  return(prunedSubgraphs)
}

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
