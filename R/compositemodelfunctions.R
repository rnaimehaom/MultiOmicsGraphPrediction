#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param pairs A list of pairs to include in the composite model.
#' @param modelResults A ModelResults object.
#' @param verbose Whether or not to print out each step.
#' @param makePlots Whether or not to plot the pruned model at each step.
#' @param pruningMethod Set to "information.gain", "odds.ratio", or "error.t.test"
#' @param binCount The number of bins to divide your original data into. Default is 10.
#' Used in information gain pruning only.
#' @param margin The margin of error for a prediction to be considered "correct".
#' Default is 0.1. Used in odds ratio pruning only.
#' @return A final predicted value for each sample in the input data
#' @export
DoSignificancePropagation <- function(pairs, modelResults, covar = NULL, verbose = FALSE, makePlots = FALSE,
                                      pruningMethod = "odds.ratio", binCount = 10, margin = 0.1){
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
                                                              makePlots = makePlots,
                                                              pruningMethod = pruningMethod,
                                                              binCount = binCount,
                                                              margin = margin)
    consolidated <- MultiOmicsGraphPrediction::ObtainCompositeModels(pairsInEachPredictor = consolidated$expandedCompositeModels, 
                                                                     importantPairs = prunedPairs)
    prevModels <- prunedPairs
  }
  
  # If number of composite models > 1, then concatenate all models and prune.
  if(length(consolidated$compositeModels) > 1){
    consolidated$compositeModels <- list(unlist(consolidated$compositeModels))
    consolidated$expandedCompositeModels <- list(unlist(consolidated$expandedCompositeModels))
    consolidated$mapping$to <- rep(1, nrow(consolidated$mapping))
    prunedPairs <- MultiOmicsGraphPrediction::PrunePredictors(compositeSubgraphs = consolidated,
                                                              previousModels = prevModels,
                                                              modelResults = modelResults, 
                                                              verbose = verbose,
                                                              makePlots = makePlots,
                                                              pruningMethod = pruningMethod,
                                                              binCount = binCount,
                                                              margin = margin)
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
  mask <- modelResults@model.input@mask
  if(nrow(weights) == 1){mask <- t(mask)}
  
  # Scale weights to sum to 1.
  #if(nrow(weights) > 1 && length(pairs) > 1){
  #  weights[,pairs] <- weights[,pairs] / matrix(rep(rowSums(weights[,pairs]), length(pairs)), ncol = length(pairs))
  #}else{
  #  weights[,pairs] <- weights[,pairs] / sum(weights[,pairs])
  #}
  

  # Analyte 1
  weighted_a1 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_a1 <- weighted_a1 + weights[,pair] * mask[,pair] * analyte1Vals[strsplit(pair, "__")[[1]][1],]
  }
  
  # beta0
  weighted_sum_b0 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b0 <- weighted_sum_b0 + weights[,pair] * mask[,pair] * rep(covariates[pair,"(Intercept)"], nrow(weights))
  }
  
  # beta1
  weighted_sum_b1 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b1 <- weighted_sum_b1 + weights[,pair] * mask[,pair] * rep(covariates[pair, "a"], nrow(weights)) * 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # beta2
  weighted_sum_b2 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b2 <- weighted_sum_b2 + weights[,pair] * mask[,pair] * rep(covariates[pair, "type"], nrow(weights))
  }
  
  # beta3
  weighted_sum_b3 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b3 <- weighted_sum_b3 + weights[,pair] * mask[,pair] * rep(covariates[pair, "a:type"], nrow(weights))* 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # covariates
  weighted_sum_covars <- rep(0, nrow(weights))
  if(!is.null(covar)){
    weighted_sum_each <- lapply(covar, function(c){
      weighted_sum <- rep(0, nrow(weights))
      for(i in 1:length(pairs)){
        pair <- pairs[[i]]
        weighted_sum <- weighted_sum + weights[,pair] * mask[,pair] * rep(covariates[pair, c], nrow(weights)) *
          covariateVals[,c]
      }
      return(weighted_sum)
    })
    weighted_sum_covars <- Reduce('+', weighted_sum_each)
  }
  
  # Analyte 1
  weighted_a1_nomask <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_a1_nomask <- weighted_a1_nomask + weights[,pair] * analyte1Vals[strsplit(pair, "__")[[1]][1],]
  }
  
  # beta0
  weighted_sum_b0_nomask <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b0_nomask <- weighted_sum_b0_nomask + weights[,pair] * rep(covariates[pair,"(Intercept)"], nrow(weights))
  }
  
  # beta1
  weighted_sum_b1_nomask <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b1_nomask <- weighted_sum_b1_nomask + weights[,pair] * rep(covariates[pair, "a"], nrow(weights)) * 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # beta2
  weighted_sum_b2_nomask <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b2_nomask <- weighted_sum_b2_nomask + weights[,pair] * rep(covariates[pair, "type"], nrow(weights))
  }
  
  # beta3
  weighted_sum_b3_nomask <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_sum_b3_nomask <- weighted_sum_b3_nomask + weights[,pair] * rep(covariates[pair, "a:type"], nrow(weights))* 
      analyte2Vals[strsplit(pair, "__")[[1]][2],]
  }
  
  # covariates
  weighted_sum_covars_nomask <- rep(0, nrow(weights))
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
    weighted_sum_covars_nomask <- Reduce('+', weighted_sum_each)
  }
  
  # Final value.
  denom <- weighted_sum_b2 + weighted_sum_b3
  denom[which(denom == 0)] <- 0.0001
  final_val <- (weighted_a1 - weighted_sum_b0 - weighted_sum_b1 - weighted_sum_covars) / denom
  
  denom_nomask <- weighted_sum_b2_nomask + weighted_sum_b3_nomask
  denom_nomask[which(denom_nomask == 0)] <- 0.0001
  final_val_nomask <- (weighted_a1_nomask - weighted_sum_b0_nomask - weighted_sum_b1_nomask - weighted_sum_covars_nomask) / denom_nomask
  #print(paste(weighted_a1_nomask - weighted_sum_b0_nomask - weighted_sum_b1_nomask - weighted_sum_covars_nomask, denom_nomask))
  #print(sum(weights[,pairs]))
  #print(paste(weighted_sum_b2_nomask[which(final_val_nomask < -5)], weighted_sum_b3_nomask[which(final_val_nomask < -5)],
  #            (weighted_a1 - weighted_sum_b0 - weighted_sum_b1 - weighted_sum_covars)[which(final_val_nomask < -5)]))
  # print(paste(weighted_sum_b2_nomask[which(final_val_nomask > 5)], weighted_sum_b3_nomask[which(final_val_nomask > 5)],
  #             (weighted_a1 - weighted_sum_b0 - weighted_sum_b1 - weighted_sum_covars)[which(final_val_nomask > 5)]))
  #print(paste(weighted_sum_b2[which(final_val > 5)], weighted_sum_b3[which(final_val > 5)],
  #            (weighted_a1 - weighted_sum_b0 - weighted_sum_b1 - weighted_sum_covars)[which(final_val > 5)]))
  #print(final_val)
  #final_val_nomask[which(final_val_nomask < -5)] <- -5
  #final_val_nomask[which(final_val_nomask > 5)] <- 5
  return(final_val_nomask)
}

#' Compute the significance value for a given prediction. You may use information
#' gain, odds ratio, or t-statistic.
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param pruningMethod Set to "information.gain", "odds.ratio", or "error.t.test"
#' @param binCount The number of bins to divide your original data into. Default is 10.
#' @param margin The margin of error for a prediction to be considered "correct".
#' Default is 0.1
#' @export
ComputeSignificance <- function(pred, trueVal, pruningMethod = "odds.ratio",
                                binCount = 10, margin = 0.1){
  if(pruningMethod == "information.gain"){
    return(ComputeInfoGain(pred = pred, trueVal = trueVal, binCount = binCount))
  }else if(pruningMethod == "odds.ratio"){
    return(ComputeOddsRatio(pred = pred, trueVal = trueVal, margin = margin))
  }else if(pruningMethod == "error.t.test"){
    return(ComputeTScore(pred = pred, trueVal = trueVal))
  }else{
    stop(paste(pruningMethod, "is not a valid pruning method."))
  }
}

#' Compute t-score for a prediction using the paired t-test. The t-score
#' tells us how far away the prediction errors are from the mean errors.
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @export
ComputeTScore <- function(pred, trueVal){
  trueVal <- trueVal[which(!is.nan(pred))]
  pred <- pred[which(!is.nan(pred))]
  tScore <- NA
  
  if(length(trueVal)>0 && length(pred)>0){
    
    # Compute the absolute error values.
    meanTrue <- rep(mean(trueVal), length(pred))
    meanAbsError <- abs(trueVal - meanTrue)
    predAbsError <- abs(trueVal - pred)
    
    # Find the difference between mean errors and prediction errors.
    # https://www.cuemath.com/data/paired-t-test/
    errorDiff <- meanAbsError - predAbsError
    n <- length(errorDiff)
    nSse <- n * sum(errorDiff ^ 2)
    errorDiffSquared <- sum(errorDiff) ^ 2
    tScore <- sum(errorDiff) / sqrt((nSse - errorDiffSquared) / (n-1))
  }
  # If all terms are NA, then categorize this as a very low t-score.
  else{
    tScore <- -1000
  }
  print(sum(errorDiff))
  print(sqrt((nSse - errorDiffSquared) / (n-1)))
  return(tScore)
}

#' Compute the odds for a given prediction. The odds ratio
#' is computed as the odds of predicting correctly (within a margin) using
#' the subgraph against correct prediction using the mean true prediction.
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param margin The margin of error for a prediction to be considered "correct".
#' Default is 0.1
#' @export
ComputeOddsRatio <- function(pred, trueVal, margin = 0.1){
  trueVal <- trueVal[which(!is.nan(pred))]
  pred <- pred[which(!is.nan(pred))]
  odds <- NA

  if(length(trueVal)>0 && length(pred)>0){
    
    # Compute the odds that the mean prediction is within the margin of error.
    meanTrue <- rep(mean(trueVal), length(pred))
    meanAbsError <- abs(trueVal - meanTrue)
    oddsForMean <- length(which(meanAbsError < margin)) / length(meanAbsError)
    
    # Compute the odds that the prediction is within the margin of error.
    predAbsError <- abs(trueVal - pred)
    oddsForPred <- length(which(predAbsError < margin)) / length(predAbsError)
    
    # Compute the odds ratio. 
    odds <- oddsForPred / oddsForMean
  }
  # If all terms are NA, then categorize this as 0 odds of correct prediction.
  else{
    odds <- 0
  }
  return(odds)
}

#' Compute the information gain for a given prediction. The information gain
#' is computed as the difference between the entropy of the true values and the
#' conditional entropy of the true value given the prediction. 
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param binCount The number of bins to divide your original data into. Default is 10.
#' @export
ComputeInfoGain <- function(pred, trueVal, binCount = 10){
  trueVal <- trueVal[which(!is.nan(pred))]
  pred <- pred[which(!is.nan(pred))]
  IG <- NULL
  
  if(length(trueVal)>0 && length(pred)>0){
    
    # Get bin size of true values given bin count. Split into bins.
    binSize <- (max(trueVals) - min(trueVals)) / binCount
    trueBins <- seq(min(trueVals), max(trueVals) + binSize, binSize)
    
    # Compute probability of true value taking each possible bin value.
    trueProbs <- hist(trueVals, breaks = trueBins, plot = FALSE)$density
    trueProbs <- trueProbs / sum(trueProbs)
    
    # Compute original entropy. Zero probability terms will be NA.
    entropyTerms <- unlist(lapply(trueProbs, function(p){
      return(-1 * p * log2(p))
    }))
    originalEntropy <- sum(entropyTerms, na.rm = TRUE)
    
    # Split the predictions into bins and find the probability of the prediction
    # taking each possible bin value.
    predBins <- seq(min(pred), max(pred) + binSize, binSize)
    predProbs <- hist(pred, breaks = predBins, plot = FALSE)$density
    predProbs <- predProbs / sum(predProbs)

    # Find the joint probability of each true and predicted bin value combination.
    jointProbs <- as.matrix(do.call(cbind, lapply(1:(length(trueBins)-1), function(tb){
      
      # Find the joint probability with each predicted bin, given the true bin.
      whichInBin <- intersect(which(trueVals >= trueBins[tb]), which(trueVals < trueBins[tb + 1]))
      predProbPerBin <- unlist(lapply(1:(length(predBins)-1), function(pb){
        whichInProbBin <- intersect(which(pred >= predBins[pb]), which(pred < predBins[pb + 1]))
        probInBothBins <- length(intersect(whichInBin, whichInProbBin)) / length(pred)
        return(probInBothBins)
      }))
      return(as.data.frame(predProbPerBin))
    })))
    
    # Compute joint entropy.
    predProbVec <- matrix(rep(predProbs, length(trueProbs)), ncol = length(trueProbs))
    # Due to rounding errors, we may occasionally get values > 0. This is fixed by clipping at 0.
    logTerm <- pmin(0, log2(jointProbs / predProbVec))
    entropyTerms <- -1 * jointProbs * logTerm
    positionsOfInterest <- intersect(which(!is.nan(entropyTerms)), which(!is.infinite(entropyTerms)))
    entropyTermsToAdd <- entropyTerms[positionsOfInterest]
    conditionalEntropy <- sum(entropyTermsToAdd, na.rm = TRUE)

    # Compute information gain.
    IG <- originalEntropy - conditionalEntropy
  }
  # If all terms are NA, then categorize this as an information loss.
  else{
    IG <- -1
  }
  return(IG)
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
#' @param pruningMethod Set to "information.gain", "odds.ratio", or "error.t.test"
#' @param binCount The number of bins to divide your original data into. Default is 10.
#' Used in information gain pruning only.
#' @param margin The margin of error for a prediction to be considered "correct".
#' Default is 0.1. Used in odds ratio pruning only.
#' @return A list of informative pairs
#' @export
PrunePredictors <- function(compositeSubgraphs, previousModels, modelResults, verbose = FALSE, makePlots = FALSE,
                            pruningMethod = "odds.ratio", binCount = 10, margin = 0.1){
  # Extract relevant information.
  pairs <- compositeSubgraphs$compositeModels
  mapping <- compositeSubgraphs$mapping
  
  prunedSubgraphs <- lapply(1:length(pairs), function(i){
    # Print initial graph.
    if(verbose == TRUE){
      print(paste("subgraph", i, ":", paste(sort(unlist(pairs[[i]])), collapse = ", ")))
    }
    
    # Initialize the predictor to include all pairs in the composite subgraph.
    importantPairs <- pairs[[i]]
    compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = importantPairs, 
                                                                     modelResults = modelResults)
    significance <- MultiOmicsGraphPrediction::ComputeSignificance(pred = unlist(compositeModel), 
                                                           trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
                                                           pruningMethod = pruningMethod,
                                                           binCount = binCount,
                                                           margin = margin)
    if(verbose == TRUE){
      print(paste(list("Original", pruningMethod, "is", significance), collapse = " "))
    }

    removedLastTime <- FALSE
    compositeModelFull <- compositeModel
    significanceFull <- significance
    
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
        significanceFull <- significance
      }
      
      # Compute the significance for the new model.
      modelsToInclude <- setdiff(importantModels, model)
      pairsToInclude <- setdiff(importantPairs, modelPairs)
      if(length(pairsToInclude) > 0){
        compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = pairsToInclude,
                                                                         modelResults = modelResults)
        significance <- MultiOmicsGraphPrediction::ComputeSignificance(pred = unlist(compositeModel),
                                                               trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
                                                               pruningMethod = pruningMethod,
                                                               binCount = binCount,
                                                               margin = margin)
      }
      
      # If the significance of the new model is greater than or equal to the full model, remove the pair.
      # If only one analyte pair remains and it has not improved over the full model, keep the pair.
      # If only one analyte pair remains and the value is 0, remove the pair.
      # If the significance of the new model is less than the full model, remove the pair.
      removedLastTime <- FALSE
      if(significance >= significanceFull && length(pairsToInclude) > 0){
        # Print the increase in significance.
        if(verbose == TRUE){
          print(paste(list("Removed model.", pruningMethod, "after removing", model, "is", significance, ", equal or greater than",
                           significanceFull), collapse = " "))
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
      }else if(significance <= 0 && length(pairsToInclude) == 0){
        if(verbose == TRUE){
          print(paste(list("Removed", model, "because final", pruningMethod, "is", significance), collapse = " "))
        }
        importantPairs <- pairsToInclude
        importantModels <- modelsToInclude
        removedLastTime <- TRUE
      }else if (length(pairsToInclude) == 0){
        if(verbose == TRUE){
          print(paste(list("Kept", model, "because it is the last remaining model"), collapse = " "))
        }
      }else{
        # Print the decrease in information gain.
        if(verbose == TRUE){
          print(paste(list("Kept model.", pruningMethod, "after removing", model, "is", significance, ", which is lower than",
                           significanceFull), collapse = " "))
        }
      }
    }
    return(importantPairs)
  })
  lengths <- unlist(lapply(prunedSubgraphs, function(g){return(length(g))}))
  prunedSubgraphs <- prunedSubgraphs[which(lengths > 0)]
  return(prunedSubgraphs)
}

#' Compute the weight of each predictor given the weights of different
#' importance metrics.
#' @param modelResults A ModelResults object.
#' @return A weight matrix for each sample and each predictor.
#' @export
ComputeImportanceWeights <- function(modelResults){
  weights <- lapply(1:length(modelResults@model.input@importance), function(i){
    imp <- modelResults@model.input@importance[[i]] * modelResults@current.importance.weights[i]
    return(as.matrix(imp))
  })
  weights_all <- Reduce("+",weights)
  #sum_weights_all <- matrix(rep(rowSums(weights_all), ncol(weights_all)), ncol = ncol(weights_all))
  #weights_all <- weights_all / sum_weights_all
  return(weights_all)
}
