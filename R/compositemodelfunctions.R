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
#' @param includeVarianceTest Scale the t-score by the f-statistic (the ratio of variances).
#' Only applicable when the pruning method is error.t.test. Default is FALSE.
#' @param modelRetention Strategy for model retention. "stringent" (the default)
#' retains only models that improve the prediction score. "lenient" also retains models that
#' neither improve nor reduce the prediction score.
#' @return A final predicted value for each sample in the input data
#' @export
DoSignificancePropagation <- function(pairs, modelResults, covar = NULL, verbose = FALSE, makePlots = FALSE,
                                      pruningMethod = "odds.ratio", binCount = 10, margin = 0.1,
                                      includeVarianceTest = FALSE, modelRetention = "stringent"){
  # Initialize consolidated pairs.
  consolidated <- list(compositeModels = pairs, expandedCompositeModels = pairs, 
                       mapping = data.frame(from = 1:length(unlist(pairs)),
                                            to = unlist(lapply(1:length(pairs), function(i){rep(i, length(pairs[[i]]))}))))
  prevModels <- unlist(pairs)
  
  # If this is the first iteration or if the mapping has not changed, stop.
  while(length(unique(consolidated$mapping$from)) != length(unique(consolidated$mapping$to))){
    prunedPairs <- MultiOmicsGraphPrediction::PrunePredictors(compositeSubgraphs = consolidated,
                                                              previousModels = prevModels,
                                                              modelResults = modelResults, 
                                                              verbose = verbose,
                                                              makePlots = makePlots,
                                                              pruningMethod = pruningMethod,
                                                              binCount = binCount,
                                                              margin = margin,
                                                              includeVarianceTest = includeVarianceTest,
                                                              modelRetention = modelRetention)

    consolidated <- MultiOmicsGraphPrediction::ObtainCompositeModels(pairsInEachPredictor = consolidated$expandedCompositeModels, 
                                                                     importantModels = prunedPairs)
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
                                                              margin = margin,
                                                              includeVarianceTest = includeVarianceTest,
                                                              modelRetention = modelRetention)
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
  analyteTgtVals <- modelResults@model.input@input.data@analyteType1
  if(modelResults@model.input@outcome == 2){
    analyteTgtVals <- modelResults@model.input@input.data@analyteType2
  }
  analyteSrcVals <- modelResults@model.input@input.data@analyteType2
  if(modelResults@model.input@independent.var.type == 1){
    analyteSrcVals <- modelResults@model.input@input.data@analyteType1
  }
  covariateVals <- modelResults@model.input@input.data@sampleMetaData
  weights <- ComputeMetaFeatureWeights(modelResults = modelResults,
                                       metaFeatures = modelResults@model.input@metaFeatures)

  # Analyte 1
  weighted_a1 <- rep(0, nrow(weights))
  for(i in 1:length(pairs)){
    pair <- pairs[[i]]
    weighted_a1 <- weighted_a1 + weights[,pair] * analyteTgtVals[strsplit(pair, "__")[[1]][2],]
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
      analyteSrcVals[strsplit(pair, "__")[[1]][1],]
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
      analyteSrcVals[strsplit(pair, "__")[[1]][1],]
  }
  
  # covariates
  weighted_sum_covars <- rep(0, nrow(weights))
  if(length(covar) > 0){
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
  
  denom <- weighted_sum_b2 + weighted_sum_b3
  denom[which(denom == 0)] <- 0.0001
  final_val <- (weighted_a1 - weighted_sum_b0 - weighted_sum_b1 - weighted_sum_covars) / denom
  return(final_val)
}

#' Compute the significance value for a given prediction. You may use information
#' gain, odds ratio, or t-statistic.
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param pruningMethod Set to "information.gain", "odds.ratio", or "error.t.test"
#' @param binCount The number of bins to divide your original data into. Default is 10.
#' @param margin The margin of error for a prediction to be considered "correct".
#' Default is 0.1
#' @param includeVarianceTest Scale the t-score by the f-statistic (the ratio of variances).
#' Only applicable when the pruning method is error.t.test. Default is FALSE.
#' @export
ComputeSignificance <- function(pred, trueVal, pruningMethod = "odds.ratio",
                                binCount = 10, margin = 0.1, includeVarianceTest = FALSE){
  if(pruningMethod == "information.gain"){
    return(ComputeInfoGain(pred = pred, trueVal = trueVal, binCount = binCount))
  }else if(pruningMethod == "odds.ratio"){
    return(ComputeOddsRatio(pred = pred, trueVal = trueVal, margin = margin))
  }else if(pruningMethod == "error.t.test"){
    return(ComputeTScore(pred = pred, trueVal = trueVal, includeVarianceTest = includeVarianceTest))
  }else{
    stop(paste(pruningMethod, "is not a valid pruning method."))
  }
}

#' Compute t-score for a prediction using the paired t-test. The t-score
#' tells us how far away the prediction errors are from the mean errors.
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param includeVarianceTest Scale the t-score by the f-statistic (the ratio of variances).
#' Default is FALSE.
#' @export
ComputeTScore <- function(pred, trueVal, includeVarianceTest = FALSE){
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
    # errorDiff <- meanAbsError - predAbsError
    # n <- length(errorDiff)
    # nSse <- n * sum(errorDiff ^ 2)
    # errorDiffSquared <- sum(errorDiff) ^ 2
    # tScore <- sum(errorDiff) / sqrt((nSse - errorDiffSquared) / (n-1))
    
    # Compute Welch's t-test.
    n <- length(meanAbsError)
    meanMeanAbsError <- mean(meanAbsError)
    meanPredAbsError <- mean(predAbsError)
    stdevMeanAbsError <- sd(meanAbsError) / sqrt(n)
    stdevPredAbsError <- sd(predAbsError) / sqrt(n)
    tScore <- (meanMeanAbsError - meanPredAbsError) / sqrt((stdevMeanAbsError ^ 2) + (stdevPredAbsError ^ 2))
    
    # If variance test is also included, compute f-statistic and add.
    # if(includeVarianceTest == TRUE){
    #   fScore <- var(predAbsError) / var(meanAbsError)
    #   tScoreOld <- tScore
    #   tScore <- tScore - fScore
    #   if(tScore > 7.93){
    #     print(paste(tScore, sum(errorDiff), sqrt((nSse - errorDiffSquared) / (n-1))))
    #   }
    # }
  }
  # If all terms are NA, then categorize this as a very low t-score.
  else{
    tScore <- -1000
  }
  
  # Plot histogram.
  # predHist <- hist(predAbsError, plot = FALSE, breaks = 10)
  # meanHist <- hist(meanAbsError, plot = FALSE, breaks = 10)
  # overallMin <- pmin(min(predAbsError), min(meanAbsError))
  # overallMax <- pmax(max(predAbsError), max(predAbsError))
  # overallMaxHeight <- length(predAbsError)
  # plot(meanHist, col=rgb(0,0,1,1/4), xlim=c(overallMin, overallMax), ylim = c(0, overallMaxHeight),
  #      main = tScore)
  # plot(predHist, col=rgb(1,0,0,1/4), add=T)
  
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

#' Obtain the list of pairs in a composite model. Do this by merging the two
#' composite models with the most overlap. Use the full models so that we
#' consider the original overlap before pruning. This uses only the structure
#' of the graph and requires no dependencies on previous layers.
#' @param pairsInEachPredictor A list of pairs contained within each predictor.
#' @param importantModels A list of all pairs that were found to be important in the
#' previous layer.
#' @export
ObtainCompositeModels <- function(pairsInEachPredictor, importantModels){
  
  # Initialize.
  compositePredictors <- lapply(importantModels, function(model){
    return(unlist(model))
  })
  compositeFull <- pairsInEachPredictor
  pairMapping <- data.frame(from = 1:length(importantModels), to = 1:length(importantModels))
  
  # Only merge composite models if there is more than one composite model to begin with.
  if(length(importantModels) > 1){
    # Compute a matrix of composite model overlaps.
    overlapsList <- lapply(1:(length(pairsInEachPredictor)-1), function(i){
      
      # Initialize model.
      model1 <- pairsInEachPredictor[[i]]
      
      # Calculate the overlap of this model with all other models.
      overlapsWithModel1 <- unlist(lapply(1:length(pairsInEachPredictor), function(j){
        model2 <- pairsInEachPredictor[[j]]
        overlap <- NA
        if(i < j){
          overlap <- length(intersect(model1, model2)) / pmin(length(model1), length(model2))
        }
        return(overlap)
      }))
      
      # Convert to data frame and return.
      overlapsWithModel1df <- as.data.frame(overlapsWithModel1)
      rownames(overlapsWithModel1df) <- 1:length(pairsInEachPredictor)
      colnames(overlapsWithModel1df) <- i
      return(overlapsWithModel1df)
    })
    overlaps <- do.call(cbind, overlapsList)
    
    # Find maximum overlap pair.
    unlistedOverlaps <- unlist(overlaps)
    combs <- expand.grid(rownames(overlaps), colnames(overlaps))
    names(unlistedOverlaps) <- paste(combs$Var2, combs$Var1, sep = ".")
    unlistedOverlaps <- unlistedOverlaps[which(!is.na(unlistedOverlaps))]
    pairToMerge <- as.numeric(strsplit(names(unlistedOverlaps)[which.max(unlistedOverlaps)], split = ".", fixed = T)[[1]])
    
    # Replace old predictor information with new.
    newPredictor <- unique(unlist(compositePredictors[pairToMerge]))
    newFull <- unique(unlist(compositeFull[pairToMerge]))
    compositePredictors[[min(pairToMerge)]] <- newPredictor
    compositeFull[[min(pairToMerge)]] <- newFull
    compositePredictors <- compositePredictors[-max(pairToMerge)]
    compositeFull <- compositeFull[-max(pairToMerge)]

    # Update the mappings.
    pairMapping[which(pairMapping$from %in% pairToMerge), "to"] <- min(pairToMerge)

    # Reassign mapping numbers.
    seq <- pairMapping$to
    compositePredictorIds <- sort(unique(pairMapping$to))
    pairMappingCopy <- pairMapping
    for(num in 1:length(compositePredictorIds)){
      pairMapping$to[which(pairMappingCopy$to == compositePredictorIds[num])] <- num
    }
  }
  
  # Return the data.
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
#' @param includeVarianceTest Scale the t-score by the f-statistic (the ratio of variances).
#' Only applicable when the pruning method is error.t.test. Default is FALSE.
#' @param tolerance Tolerance factor when computing equality of two numeric values.
#' @param modelRetention Strategy for model retention. "stringent" (the default)
#' retains only models that improve the prediction score. "lenient" also retains models that
#' neither improve nor reduce the prediction score.
#' @return A list of informative pairs
#' @export
PrunePredictors <- function(compositeSubgraphs, previousModels, modelResults, verbose = FALSE, makePlots = FALSE,
                            pruningMethod = "odds.ratio", binCount = 10, margin = 0.1,
                            includeVarianceTest = includeVarianceTest, tolerance = 1e-5,
                            modelRetention = "stringent"){
  
  # Extract relevant information.
  pairs <- compositeSubgraphs$compositeModels
  mapping <- compositeSubgraphs$mapping
  
  prunedSubgraphs <- lapply(1:length(pairs), function(i){
    # Print initial graph.
    if(verbose == TRUE){
      print(paste("subgraph", i, ":", paste(sort(unlist(pairs[[i]])), collapse = ", ")))
    }
    
    # Initialize the predictor to include all pairs in the composite subgraph.
    compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = pairs[[i]], 
                                                                     modelResults = modelResults)
    significance <- MultiOmicsGraphPrediction::ComputeSignificance(pred = unlist(compositeModel), 
                                                                   trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
                                                                   pruningMethod = pruningMethod,
                                                                   binCount = binCount,
                                                                   margin = margin,
                                                                   includeVarianceTest = includeVarianceTest)
    if(verbose == TRUE){
      print(paste(list("Original", pruningMethod, "is", significance), collapse = " "))
    }
    
    removedLastTime <- FALSE
    compositeModelFull <- compositeModel
    significanceFull <- significance
    
    # Figure out which of the previous models mapped to this one.
    previousModelsMapped <- mapping$from[which(mapping$to == i)]
    importantModels <- previousModelsMapped
    importantPairs <- previousModels[previousModelsMapped]

    # Sort the models in order of their individual performance.
    individualPerformance <- unlist(lapply(1:length(previousModelsMapped), function(m){
      model <- importantModels[m]
      modelPairs <- unlist(importantPairs[m])
      
      compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = modelPairs,
                                                                       modelResults = modelResults)
      significance <- MultiOmicsGraphPrediction::ComputeSignificance(pred = unlist(compositeModel),
                                                                     trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
                                                                     pruningMethod = pruningMethod,
                                                                     binCount = binCount,
                                                                     margin = margin,
                                                                     includeVarianceTest = includeVarianceTest)
      return(significance)
    }))
    modelRemovalOrder <- order(individualPerformance)
    
    # Sequentially test removal of each model.
    for(m in modelRemovalOrder){
      
      # Save the information gain for the full model.
      if(removedLastTime == TRUE){
        compositeModelFull <- compositeModel
        significanceFull <- significance
      }
      
      # Compute the significance for the new model.
      modelsToInclude <- setdiff(importantModels, previousModelsMapped[m])
      pairsToInclude <- intersect(unlist(importantPairs), unlist(previousModels[modelsToInclude]))
      if(length(pairsToInclude) > 0){
        compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = pairsToInclude,
                                                                         modelResults = modelResults)
        significance <- MultiOmicsGraphPrediction::ComputeSignificance(pred = unlist(compositeModel),
                                                                       trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
                                                                       pruningMethod = pruningMethod,
                                                                       binCount = binCount,
                                                                       margin = margin,
                                                                       includeVarianceTest = includeVarianceTest)
      }
      
      # If the significance of the new model is greater than or equal to the full model, remove the pair.
      # If only one analyte pair remains and it has not improved over the full model, keep the pair.
      # If only one analyte pair remains and the value is 0, remove the pair.
      # If the significance of the new model is less than the full model, remove the pair.
      significance <- round(significance, 2)
      significanceFull <- round(significanceFull, 2)
      model <- intersect(unlist(importantPairs), unlist(previousModels[setdiff(importantModels, modelsToInclude)]))
      removedLastTime <- FALSE
      meetsCutoffForRemoval <- significance >= significanceFull - tolerance
      if(modelRetention == "lenient"){
        meetsCutoffForRemoval <- significance > significanceFull
      }
      if(meetsCutoffForRemoval == TRUE && length(pairsToInclude) > 0){
        # Print the increase in significance.
        if(verbose == TRUE){
          print(paste(list("Removed model.", pruningMethod, "after removing", paste(model, collapse = ","), "is", 
                           format(significance, nsmall = 2), ", as compared to",
                           format(significanceFull, nsmall = 2)), collapse = " "))
        }
        
        importantPairs <- pairsToInclude
        importantModels <- modelsToInclude
        removedLastTime <- TRUE
      }else if (length(pairsToInclude) == 0){
        if(verbose == TRUE){
          print(paste(list("Kept", paste(model, collapse = ","), "because it is the last remaining model"), collapse = " "))
        }
      }else{
        # Print the decrease in information gain.
        if(verbose == TRUE){
          print(paste(list("Kept model.", pruningMethod, "after removing", paste(model, collapse = ","), "is", 
                           format(significance, nsmall = 2), ", as compared to",
                           format(significanceFull, nsmall = 2)), collapse = " "))
        }
      }
    }
    if(verbose == TRUE){
      print(paste("The pruned set of pairs is:", paste(unlist(importantPairs), collapse = ",")))
    }
    return(unlist(importantPairs))
  })
  lengths <- unlist(lapply(prunedSubgraphs, function(g){return(length(g))}))
  prunedSubgraphs <- prunedSubgraphs[which(lengths > 0)]
  return(prunedSubgraphs)
}

#' Compute the weight of each predictor given the weights of different
#' metafeatures.
#' @param modelResults A ModelResults object.
#' @param metaFeatures A set of metafeatures, such as that found within ModelResults
#' @return A weight matrix for each sample and each predictor.
#' @export
ComputeMetaFeatureWeights <- function(metaFeatures, modelResults){
  weights <- lapply(1:length(modelResults@model.input@metaFeatures), function(i){
    imp <- metaFeatures[[i]] * modelResults@current.metaFeature.weights[i]
    return(as.matrix(imp))
  })
  weights_all <- Reduce("+",weights)
  return(weights_all)
}
