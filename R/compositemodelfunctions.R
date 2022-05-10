#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param pairs A list of pairs to include in the composite model.
#' @param modelResults A ModelResults object.
#' @param verbose Whether or not to print out each step.
#' @param makePlots Whether or not to plot the pruned model at each step.
#' @return A final predicted value for each sample in the input data
#' @export
DoSignificancePropagation <- function(pairs, modelResults, covar = NULL, verbose = FALSE, makePlots = FALSE){
  # Initialize consolidated pairs.
  consolidated <- list(compositeModels = pairs, expandedCompositeModels = pairs, 
                       mapping = data.frame(from = 1:length(unlist(pairsPredAll)),
                                            to = unlist(lapply(1:length(pairsPredAll), function(i){rep(i, length(pairsPredAll[[i]]))}))))
  prevModels <- unlist(pairs)
  
  # If this is the first iteration or if the mapping has not changed, stop.
  while(length(unique(consolidated$mapping[,1])) != length(unique(consolidated$mapping[,2]))){
    prunedPairs <- MultiOmicsGraphPrediction::PrunePredictors(compositeSubgraphs = consolidated,
                                                              previousModels = prevModels,
                                                              modelResults = modelResults, 
                                                              verbose = verbose,
                                                              makePlots = makePlots)
    consolidated <- MultiOmicsGraphPrediction::ObtainCompositeModels(pairsInEachPredictor = consolidated$expandedCompositeModels, 
                                                                     importantPairs = prunedPairs)
    prevModels <- prunedPairs
  }
  
  # Return consolidated pairs.
  return(prunedPairs)
}

#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param pairs A list of pairs to include in the composite model.
#' @param IntLIMresults The filtered results obtained from the IntLIM function
#' ProcessRe
#' @param modelResults A ModelResults object.
#' @param covar Covariates to be used in the model.
#' @return A final predicted value for each sample in the input data
#' @export
CompositePrediction <- function(pairs, modelResults){
  
  # Set variables for further analysis.
  covar <- modelResults@model.input@covariates
  covariates <- modelResults@model.input@model.properties
  analyte1Vals <- modelResults@model.input@input.data@analyteType1
  analyte2Vals <- modelResults@model.input@input.data@analyteType2
  covariateVals <- modelResults@model.input@input.data@sampleMetaData
  weights <- ComputeImportanceWeights(modelResults)

  # Analyte 1
  weighted_a1 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_a1 <- weighted_a1 + weights[,pair] * analyte1Vals[strsplit(pair, "__")[[1]][1],]
  }
  
  # beta0
  weighted_sum_b0 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b0 <- weighted_sum_b0 + weights[,pair] * rep(covariates[pair,"(Intercept)"], nrow(weights))
  }
  
  # beta1
  weighted_sum_b1 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b1 <- weighted_sum_b1 + weights[,pair] * rep(covariates[pair, "a"], nrow(weights)) * 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # beta2
  weighted_sum_b2 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b2 <- weighted_sum_b2 + weights[,pair] * rep(covariates[pair, "type"], nrow(weights))
  }
  
  # beta3
  weighted_sum_b3 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b3 <- weighted_sum_b3 + weights[,pair] * rep(covariates[pair, "a:type"], nrow(weights))* 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # covariates
  weighted_sum_covars <- rep(0, nrow(weights))
  if(!is.null(covar)){
    weighted_sum_each <- lapply(covar, function(c){
      weighted_sum <- rep(0, nrow(weights))
      for(i in 1:length(pairs)){
        pair <- pairs[[i]]
        weighted_sum <- weighted_sum + weights[,pair] * rep(covariates[pair, c], nrow(weights)) *
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
#' @param percentOverlapCutoff Neighborhoods that have this percent overlap
#' with another neighborhood will be removed from consideration. Default is 100.
#' @export
ObtainSubgraphNeighborhoods <- function(modelInput, percentOverlapCutoff = 100){
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
      if(length(intersect(nodeAndNeighbors, nbrAndNeighbors)) /
         length(nodeAndNeighbors) * 100 >= percentOverlapCutoff){
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
ObtainCompositeModels <- function(pairsInEachPredictor, importantPairs){

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

#' Given multiple composite predictors, prune the predictors that are not needed.
#' @param compositeSubgraphs A list of pairs to include in the composite model.
#' @param previousModels A list of the previous models that were consolidated.
#' @param modelResults A ModelResults object.
#' @param verbose Whether or not to print out each step.
#' @param makePlots Whether or not to plot the pruned model at each step.
#' @return A list of informative pairs
#' @export
PrunePredictors <- function(compositeSubgraphs, previousModels, modelResults, verbose = FALSE, makePlots = FALSE){
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
                                                                     modelResults = modelResults)
    infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel), 
                                                           trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype])
    if(verbose == TRUE){
      print(paste(list("Original information gain is", infoGain), collapse = " "))
    }
    
    removedLastTime <- FALSE
    compositeModelFull <- compositeModel
    infoGainFull <- infoGain
    
    # Figure out which of the previous models mapped to this one.
    previousModelsMapped <- mapping$from[which(mapping$to == i)]
    importantModels <- previousModelsMapped
    previousPairsMapped <- previousModels[previousModelsMapped]
    
    # Sequentially test removal of each model.
    for(m in 1:length(previousModelsMapped)){
      model <- previousModelsMapped[m]
      modelPairs <- unlist(previousPairsMapped[m])
      
      # Save the information gain for the full model.
      if(removedLastTime == TRUE){
        compositeModelFull <- compositeModel
        infoGainFull <- infoGain
      }
      
      # Compute the information gain for the new model.
      modelsToInclude <- setdiff(importantModels, model)
      pairsToInclude <- setdiff(importantPairs, modelPairs)
      if(length(pairsToInclude) > 0){
        compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = pairsToInclude,
                                                                         modelResults = modelResults)
        infoGain <- MultiOmicsGraphPrediction::ComputeInfoGain(pred = unlist(compositeModel),
                                                               trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype])
      }
      
      # If the information gain for the new model is greater than for the full model, remove the pair.
      removedLastTime <- FALSE
      if(infoGain >= infoGainFull && length(pairsToInclude) > 0){
        # Print the increase in information gain.
        if(verbose == TRUE){
          print(paste(list("Removed model. Information gain after removing", model, "is", infoGain - infoGainFull), collapse = " "))
        }
        
        importantPairs <- pairsToInclude
        importantModels <- modelsToInclude
        removedLastTime <- TRUE
        
        tryCatch({
          if(makePlots == TRUE){
            # Plot the new line graph.
            newPredictions <- matrix(rep(compositeModel, ncol(modelResults@model.input@node.wise.prediction)),
                                       ncol = ncol(modelResults@model.input@node.wise.prediction))
            rownames(newPredictions) <- rownames(modelResults@model.input@node.wise.prediction)
            colnames(newPredictions) <- colnames(modelResults@model.input@node.wise.prediction)
            modelResultsNew <- modelResults
            modelResultsNew@model.input@node.wise.prediction <- newPredictions
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
#' @return A weight matrix for each sample and each predictor.
#' @export
ComputeImportanceWeights <- function(modelResults){
  importanceWeights <- modelResults@current.importance.weights
  weights <- lapply(1:length(modelResults@model.input@importance), function(i){
    imp <- modelResults@model.input@importance[[i]] * modelResults@current.importance.weights[i]
    return(as.matrix(imp))
  })
  weights_all <- Reduce("+",weights)
  return(weights_all)
}
