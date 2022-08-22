#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param pairs A list of pairs to include in the composite model.
#' @param modelResults A ModelResults object.
#' @param verbose Whether or not to print out each step.
#' @param pruningMethod The method to use for pruning. Right now, only "error.t.test" is valid.
#' @param modelRetention Strategy for model retention. "stringent" (the default)
#' retains only models that improve the prediction score. "lenient" also retains models that
#' neither improve nor reduce the prediction score.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @param weightCutoff Only consider pairs with weight above this cutoff.
#' @return A final predicted value for each sample in the input data
#' @export
DoSignificancePropagation <- function(pairs, modelResults, covar = c(), verbose = FALSE,
                                      pruningMethod = "error.t.test", modelRetention = "stringent",
                                      minCutoff, maxCutoff, useCutoff = FALSE, weightCutoff = 0){                       	
                              	
  # Compute weights for each predictor.
  weights <- ComputeMetaFeatureWeights(modelResults = modelResults,
                                       metaFeatures = modelResults@model.input@metaFeatures)

  # Initialize consolidated pairs.
  targets <- pairs$target
  sources <- pairs$source
  pairs <- pairs$edge
  consolidated <- list(compositeModels = pairs, expandedCompositeModels = pairs, 
                       compositeModelTargets = targets, expandedCompositeModelTargets = targets,
                       compositeModelSource = sources, expandedCompositeModelSources = sources,
                       mapping = data.frame(from = 1:length(unlist(pairs)),
                                            to = unlist(lapply(1:length(pairs), function(i){rep(i, length(pairs[[i]]))}))))
  prevModels <- list(pairs = unlist(pairs), targets = unlist(targets), sources = unlist(sources))

  # Compute individual performance of each predictor.
  individualPerformance <- ComputeIndividualPerformance(modelResults = modelResults,
                                                        trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
                                                        pruningMethod = pruningMethod)

  # If this is the first iteration or if the mapping has not changed, stop.
  while(length(unique(consolidated$mapping$from)) != length(unique(consolidated$mapping$to))){
    prunedPairs <- MultiOmicsGraphPrediction::PrunePredictors(compositeSubgraphs = consolidated,
                                                              previousModels = prevModels,
                                                              modelResults = modelResults, 
                                                              verbose = verbose,
                                                              pruningMethod = pruningMethod,
                                                              modelRetention = modelRetention,
                                                              minCutoff = minCutoff,
                                                              maxCutoff = maxCutoff,
                                                              useCutoff = useCutoff,
                                                              weights = weights,
                                                              individualPerformance = individualPerformance,
                                                              weightCutoff = weightCutoff)

    consolidated <- MultiOmicsGraphPrediction::ObtainCompositeModels(pairsInEachPredictor = consolidated$expandedCompositeModels,
                                                                     targetsInEachPredictor = consolidated$expandedCompositeModelTargets,
                                                                     sourcesInEachPredictor = consolidated$expandedCompositeModelSources,
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
                                                              pruningMethod = pruningMethod,
                                                              modelRetention = modelRetention,
                                                              minCutoff = minCutoff,
                                                              maxCutoff = maxCutoff,
                                                              useCutoff = useCutoff,
                                                              weights = weights,
                                                              individualPerformance = individualPerformance)
  }
  
  # Return consolidated pairs.
  return(prunedPairs)
}

#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param pairs A list of pairs to include in the composite model.
#' @param modelResults A ModelResults object.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @param weights The weight for each predictor, calculated using ComputeMetaFeatureWeights()
#' @param targets Target analytes for all pairs
#' @param sources Source analytes for all pairs
#' @return A final predicted value for each sample in the input data
#' @export
CompositePrediction <- function(pairs, modelResults, minCutoff, maxCutoff, useCutoff = FALSE, weights,
                                targets, sources){

  # Run prediction for training data.
  final_val <- Predict(pairs = pairs,
    inputData = modelResults@model.input@input.data, 
                       weights = weights, 
                       model = modelResults, 
                       minCutoff = minCutoff, 
                       maxCutoff = maxCutoff, 
                       useCutoff = useCutoff,
                      useActivation = FALSE,
    independentVarType = modelResults@model.input@independent.var.type, 
    outcomeType = modelResults@model.input@outcome,
    targets = targets,
    sources = sources)
  return(final_val)
}

#' Compute the significance value for a given prediction. You may use information
#' gain, odds ratio, or t-statistic.
#' @param pairs A list of pairs to include in the composite model.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param pruningMethod The method to use for pruning. Right now, only "error.t.test" is valid.
#' @param includeVarianceTest Scale the t-score by the f-statistic (the ratio of variances).
#' Only applicable when the pruning method is error.t.test. Default is FALSE.
#' @export
ComputeSignificance <- function(pred, trueVal, pruningMethod = "error.t.test"){
  if(pruningMethod == "error.t.test"){
    return(ComputeTScore(pred = pred, trueVal = trueVal))
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
  tScore <- -1000
  if(length(trueVal)>0 && length(pred)>0){
    
    # Compute the absolute error values.
    meanTrue <- rep(mean(trueVal), length(pred))
    meanAbsError <- abs(trueVal - meanTrue)
    predAbsError <- abs(trueVal - pred)
    
    # Compute Welch's t-test.
    n <- length(meanAbsError)
    meanMeanAbsError <- mean(meanAbsError)
    meanPredAbsError <- mean(predAbsError)
    stdevMeanAbsError <- sd(meanAbsError) / sqrt(n)
    stdevPredAbsError <- sd(predAbsError) / sqrt(n)
    tScore <- (meanMeanAbsError - meanPredAbsError) / sqrt((stdevMeanAbsError ^ 2) + (stdevPredAbsError ^ 2))
  }
  
  return(tScore)
}

#' Obtain the list of pairs in a composite model.
#' @param modelInput A ModelInput object
#' @return A list of sets, where each set is a neighborhood of nodes.
#' @export
ObtainSubgraphNeighborhoods <- function(modelInput){
  # Convert to graph.
  graph <- igraph::graph_from_adjacency_matrix(modelInput@coregulation.graph)
  
  # Return the list of edges for each node so that they do not need to be looked
  # up for each edge.
  allNodes <- igraph::as_ids(igraph::V(graph))
  print("Finding pairs containing each analyte")
  edgesPerNode <- lapply(1:length(allNodes), function(i){
    node <- allNodes[i]
    nodeTargetEdges <- c()
    nodeSourceEdges <- c()
    whichTarget <- colnames(modelInput@coregulation.graph)[which(modelInput@coregulation.graph[node,] == 1)]
    whichSource <- rownames(modelInput@coregulation.graph)[which(modelInput@coregulation.graph[,node] == 1)]
    if(length(whichTarget) > 0){
      nodeTargetEdges <- paste(node, whichTarget, sep = "|")
    }
    if(length(whichSource) > 0){
      nodeSourceEdges <- paste(whichSource, node, sep = "|")
    }
    return(union(nodeTargetEdges, nodeSourceEdges))
  })
  names(edgesPerNode) <- allNodes
  
  # For each edge, find its nodes' neighbors. Ensure that the neighborhood is not a 
  # subset of any of the neighbors' neighborhoods. This prevents subset graphs
  # from being returned.
  print("Segmenting graph into neighborhoods")
  edges <- igraph::as_ids(igraph::E(graph))
  edgesToCheckQueue <- edges
  edgeNeighborhoods <- list()
  targetNeighborhoods <- list()
  sourceNeighborhoods <- list()
  while(length(edgesToCheckQueue) > 0){
    # Select the next edge to check. Find its neighborhoods.
    edge <- edgesToCheckQueue[[1]]
    nodes <- igraph::ends(graph, edge)
    edgeAndNeighbors <- union(edgesPerNode[[nodes[1]]], edgesPerNode[[nodes[2]]])
    
    # Remove the current edge from the queue
    edgesToCheckQueue <- setdiff(edgesToCheckQueue, edgeAndNeighbors)
    
    # Add the neighborhood.
    edgeAndNeighborsFormatted <- gsub("|", "__", edgeAndNeighbors, fixed = TRUE)
    edgeNeighborhoods[[length(edgeNeighborhoods) + 1]] <- edgeAndNeighborsFormatted
    targetNeighborhoods[[length(targetNeighborhoods) + 1]] <- unname(unlist(data.frame(strsplit(edgeAndNeighborsFormatted, "__"))[2,]))
    sourceNeighborhoods[[length(sourceNeighborhoods) + 1]] <- unname(unlist(data.frame(strsplit(edgeAndNeighborsFormatted, "__"))[1,]))
  }
  
  # Return neighborhoods.
  return(list(edge = edgeNeighborhoods, target = targetNeighborhoods, source = sourceNeighborhoods))
}

#' Obtain the list of pairs in a composite model. Do this by merging all composite
#' models with overlap. Use the full models so that we
#' consider the original overlap before pruning. This uses only the structure
#' of the graph and requires no dependencies on previous layers.
#' @param graph The full graph of predictors.
#' @param pairsInEachPredictor A list of pairs contained within each predictor.
#' @param sourcesInEachPredictor A list of sources for pairs contained within each predictor.
#' @param sourcesInEachPredictor A list of targets for pairs contained within each predictor.
#' @param importantModels A list of all pairs that were found to be important in the
#' previous layer.
#' @return A list with the following elements: A list of sets comprising the trimmed
#' models, a list of sets comprising the untrimmed models, and a mapping from
#' the predictors in the previous stage to the current stage.
#' @export
ObtainCompositeModels <- function(graph, pairsInEachPredictor, targetsInEachPredictor, sourcesInEachPredictor, importantModels){
  
  # Initialize.
  compositePredictors <- lapply(importantModels$pairs, function(model){
    return(unlist(model))
  })
  compositePredictorsTargets <- lapply(importantModels$targets, function(model){
    return(unlist(model))
  })
  compositePredictorsSources <- lapply(importantModels$sources, function(model){
    return(unlist(model))
  })
  compositeFull <- pairsInEachPredictor
  compositeFullTargets <- targetsInEachPredictor
  compositeFullSources <- sourcesInEachPredictor
  pairMapping <- data.frame(from = 1:length(importantModels$pairs), to = 1:length(importantModels$pairs))
  
  # Only merge composite models if there is more than one composite model to begin with.
  if(length(importantModels$pairs) > 1){
    compositePredictors <- list()
    compositeFull <- list()
    compositePredictorTargets <- list() 
    compositeFullTargets <- list()
    compositePredictorSources <- list()
    compositeFullSources <- list() 
    
    # Find all connected components of the graph. This is identical to merging all
    # clusters that have overlap.
    components <- igraph::component_distribution(graph, mode = "strong")
    for(i in 1:length(components)){
      component <- components[[i]]
      
      # Split component into sources and targets.
      compositeFullTargets <- unname(unlist(data.frame(strsplit(component, "__"))[2,]))
      compositeFullSources <- unname(unlist(data.frame(strsplit(component, "__"))[1,]))
      
      # Find what maps to the component and update mapping.
      doesMap <- unlist(lapply(compositeFull, function(model){
        return(setdiff(model, component) == 0)
      }))
      pairMapping$to[which(doesMap)] <- i
      
      # Update composite predictors.
      whichImportant <- which(components[[i]] %in% unlist(importantModels))
      compositePredictors[length(compositePredictors) + 1] <- component[whichImportant]
      compositeFull[length(compositePredictors) + 1] <- component
      compositePredictorTargets <- compositeFullTargets[whichImportant]
      compositePredictorSources <- compositeFullSources[whichImportant]
    }
  }
    # # Compute a matrix of composite model overlaps.
    # overlapsList <- lapply(1:(length(pairsInEachPredictor)-1), function(i){
    #   
    #   # Initialize model.
    #   model1 <- pairsInEachPredictor[[i]]
    #   
    #   # Calculate the overlap of this model with all other models.
    #   overlapsWithModel1 <- unlist(lapply(1:length(pairsInEachPredictor), function(j){
    #     model2 <- pairsInEachPredictor[[j]]
    #     overlap <- 0
    #     if(i < j){
    #       overlap <- length(intersect(model1, model2)) / pmin(length(model1), length(model2))
    #     }
    #     return(overlap)
    #   }))
    #   
    #   # Convert to data frame and return.
    #   overlapsWithModel1df <- as.data.frame(overlapsWithModel1)
    #   rownames(overlapsWithModel1df) <- 1:length(pairsInEachPredictor)
    #   colnames(overlapsWithModel1df) <- i
    #   return(overlapsWithModel1df)
    # })
    # overlaps <- do.call(cbind, overlapsList)

    # Find maximum overlap pair.
  #   if(max(overlaps) > 0){
  #     unlistedOverlaps <- unlist(overlaps)
  #     combs <- expand.grid(rownames(overlaps), colnames(overlaps))
  #     names(unlistedOverlaps) <- paste(combs$Var2, combs$Var1, sep = ".")
  #     unlistedOverlaps <- unlistedOverlaps[which(!is.na(unlistedOverlaps))]
  #     pairToMerge <- as.numeric(strsplit(names(unlistedOverlaps)[which.max(unlistedOverlaps)], split = ".", fixed = T)[[1]])
  #     
  #     # Replace old predictor information with new.
  #     # Figure out how to modify newPredictor and newFull for targets and sources.
  #     newPredictor <- unique(unlist(compositePredictors[pairToMerge]))
  #     newFull <- unique(unlist(compositeFull[pairToMerge]))
  #     compositePredictors[[min(pairToMerge)]] <- newPredictor
  #     compositePredictorsTargets[[min(pairToMerge)]] <- newPredictor
  #     compositePredictorsSources[[min(pairToMerge)]] <- newPredictor
  #     compositeFull[[min(pairToMerge)]] <- newFull
  #     compositeFullTargets[[min(pairToMerge)]] <- newFull
  #     compositeFullTargets[[min(pairToMerge)]] <- newFull
  #     compositePredictors <- compositePredictors[-max(pairToMerge)]
  #     compositePredictorsTargets <- compositePredictorsTargets[-max(pairToMerge)]
  #     compositePredictorsSources <- compositePredictorsSources[-max(pairToMerge)]
  #     compositeFull <- compositeFull[-max(pairToMerge)]
  #     compositeFullTargets <- compositeFullTargets[-max(pairToMerge)]
  #     compositeFullSources <- compositeFullSources[-max(pairToMerge)]
  #     
  #     # Update the mappings.
  #     pairMapping[which(pairMapping$from %in% pairToMerge), "to"] <- min(pairToMerge)
  #     
  #     # Reassign mapping numbers.
  #     seq <- pairMapping$to
  #     compositePredictorIds <- sort(unique(pairMapping$to))
  #     pairMappingCopy <- pairMapping
  #     for(num in 1:length(compositePredictorIds)){
  #       pairMapping$to[which(pairMappingCopy$to == compositePredictorIds[num])] <- num
  #     }
  #   }
  # }

  # Return the data.
  return(list(compositeModels = compositePredictors, expandedCompositeModels = compositeFull, mapping = pairMapping,
              compositeModelTargets = compositePredictorTargets, expandedCompositeModelTargets = compositeFullTargets,
              compositeModelSources = compositePredictorSources, expandedCompositeModelSources = compositeFullSources))
}

#' Given multiple composite predictors, prune the predictors that are not needed.
#' @param compositeSubgraphs A list of pairs to include in the composite model.
#' @param previousModels A list of the previous models that were consolidated.
#' @param previousModelsTarget A list of the targets in the previous models that were consolidated.
#' @param previousModelsSource A list of the sources in the previous models that were consolidated.
#' @param modelResults A ModelResults object.
#' @param verbose Whether or not to print out each step.
#' @param makePlots Whether or not to plot the pruned model at each step.
#' @param pruningMethod The method to use for pruning. Right now, only "error.t.test" is valid.
#' @param tolerance Tolerance factor when computing equality of two numeric values.
#' @param modelRetention Strategy for model retention. "stringent" (the default)
#' retains only models that improve the prediction score. "lenient" also retains models that
#' neither improve nor reduce the prediction score.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @param weights The weights for each predictor, calculated using ComputeMetaFeatureWeights()
#' @param individualPerformance The score (using the pruning method) for each individual component of the model.
#' @param weightCutoff Only consider pairs with weight above this cutoff.
#' @return A list of sets, where each set is a neighborhood of nodes.
#' @export
PrunePredictors <- function(compositeSubgraphs, previousModels, modelResults, verbose = FALSE, makePlots = FALSE,
                            pruningMethod = "error.t.test", tolerance = 1e-5,
                            modelRetention = "stringent", minCutoff, maxCutoff, useCutoff = FALSE, weights,
                            individualPerformance, weightCutoff){
  
  # Extract relevant information.
  pairs <- compositeSubgraphs$compositeModels
  targets <- compositeSubgraphs$compositeModelTargets
  sources <- compositeSubgraphs$compositeModelSource
  mapping <- compositeSubgraphs$mapping
  tstart <- Sys.time()
  prunedSubgraphs <- lapply(1:length(pairs), function(i){
    # Prune by weight.
    whichWeights <- which(apply(weights[,pairs[[i]]], 2, function(x) max(x)) > weightCutoff)
    localPairs <- pairs[[i]][whichWeights]
    localTargets <- targets[[i]][whichWeights]
    localSources <- sources[[i]][whichWeights]
    print(paste(length(localPairs), "out of", length(pairs[[i]]), "remaining after weight pruning"))
      
    # Print initial graph.
    maxToPrint = 50
    if(verbose == TRUE){
      if(length(localPairs) < maxToPrint){
        print(paste("subgraph", i, ":", paste(sort(unlist(localPairs)), collapse = ", ")))
      }else{
        print(paste("subgraph", i, ":", paste(sort(unlist(localPairs[1:maxToPrint])), collapse = ", ")))
      }
    }
    
    # Make sure there is at least one model remaining.
    if(length(localPairs) > 1){
      # Initialize the predictor to include all pairs in the composite subgraph.
      compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = localPairs, 
                                                                       targets = localTargets,
                                                                       sources = localSources,
                                                                       modelResults = modelResults,
                                                                       minCutoff = minCutoff,
                                                                       maxCutoff = maxCutoff,
                                                                       useCutoff = useCutoff,
                                                                       weights = weights)
      significance <- MultiOmicsGraphPrediction::ComputeSignificance(pred = unlist(compositeModel), 
                                                                     trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
                                                                     pruningMethod = pruningMethod)
      if(verbose == TRUE){
        print(paste(list("Original", pruningMethod, "is", significance), collapse = " "))
      }
      
    #   removedLastTime <- FALSE
    #   compositeModelFull <- compositeModel
    #   significanceFull <- significance
    #   
    #   # Figure out which of the previous models mapped to this one.
    #   previousModelsMapped <- mapping$from[which(mapping$to == i)]
    #   importantModels <- previousModelsMapped
    #   importantPairs <- previousModels$pairs[previousModelsMapped]
    #   importantTargets <- previousModels$targets[previousModelsMapped]
    #   importantSources <- previousModels$sources[previousModelsMapped]
    #   
    #   # Filter the important models.
    #   whichWeightsImportant <- which(apply(weights[,unlist(importantPairs)], 2, function(x) max(x)) > weightCutoff)
    #   importantPairs <- importantPairs[whichWeightsImportant]
    #   importantTargets <- importantTargets[whichWeightsImportant]
    #   importantSources <- importantSources[whichWeightsImportant]
    #   
    #   # Sort the models in order of their individual performance.
    #   whichPair <- which(names(individualPerformance) %in% localPairs)
    #   modelRemovalOrder <- order(individualPerformance[whichPair])
    #   
    #   # Sequentially test removal of each model.
    #   for(m in modelRemovalOrder){
    #     
    #     # Save the information gain for the full model.
    #     if(removedLastTime == TRUE){
    #       compositeModelFull <- compositeModel
    #       significanceFull <- significance
    #     }
    #     
    #     # Compute the significance for the new model.
    #     modelsToInclude <- setdiff(importantModels, previousModelsMapped[m])
    #     whichToInclude <- which(unlist(importantPairs) %in% unlist(previousModels$pairs[modelsToInclude]))
    #     pairsToInclude <- unlist(importantPairs)[whichToInclude]
    #     targetsToInclude <- unlist(importantTargets)[whichToInclude]
    #     sourcesToInclude <- unlist(importantSources)[whichToInclude]
    #     if(length(pairsToInclude) > 0){
    #       compositeModel <- MultiOmicsGraphPrediction::CompositePrediction(pairs = pairsToInclude,
    #                                                                        targets = targetsToInclude,
    #                                                                        sources = sourcesToInclude,
    #                                                                        modelResults = modelResults,
    #                                                                        minCutoff = minCutoff,
    #                                                                        maxCutoff = maxCutoff,
    #                                                                        useCutoff = useCutoff,
    #                                                                        weights = weights)
    #       significance <- MultiOmicsGraphPrediction::ComputeSignificance(pred = unlist(compositeModel),
    #                                                                      trueVal = modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype],
    #                                                                      pruningMethod = pruningMethod)
    #     }
    #     
    #     # If the significance of the new model is greater than or equal to the full model, remove the pair.
    #     # If only one analyte pair remains and it has not improved over the full model, keep the pair.
    #     # If only one analyte pair remains and the value is 0, remove the pair.
    #     # If the significance of the new model is less than the full model, remove the pair.
    #     significance <- round(significance, 2)
    #     significanceFull <- round(significanceFull, 2)
    #     model <- unlist(importantPairs)[which(unlist(importantPairs) %in% unlist(previousModels$pairs[previousModelsMapped[m]]))]
    #     removedLastTime <- FALSE
    #     meetsCutoffForRemoval <- significance >= significanceFull - tolerance
    #     if(modelRetention == "lenient"){
    #       meetsCutoffForRemoval <- significance > significanceFull
    #     }
    #     if(meetsCutoffForRemoval == TRUE && length(pairsToInclude) > 0){
    #       # Print the increase in significance.
    #       if(verbose == TRUE){
    #         print(paste(list("Removed model.", pruningMethod, "after removing", paste(model, collapse = ","), "is", 
    #                          format(significance, nsmall = 2), ", as compared to",
    #                          format(significanceFull, nsmall = 2)), collapse = " "))
    #       }
    #       
    #       importantPairs <- pairsToInclude
    #       importantTargets <- targetsToInclude
    #       importantSources <- sourcesToInclude
    #       importantModels <- modelsToInclude
    #       removedLastTime <- TRUE
    #     }else if (length(pairsToInclude) == 0){
    #       if(verbose == TRUE){
    #         print(paste(list("Kept", paste(model, collapse = ","), "because it is the last remaining model"), collapse = " "))
    #       }
    #     }else{
    #       # Print the decrease in information gain.
    #       if(verbose == TRUE){
    #         print(paste(list("Kept model.", pruningMethod, "after removing", paste(model, collapse = ","), "is", 
    #                          format(significance, nsmall = 2), ", as compared to",
    #                          format(significanceFull, nsmall = 2)), collapse = " "))
    #       }
    #     }
    #   }
      importantPairs <- localPairs
      importantTargets <- localTargets
      importantSources <- localSources
    }else{
      importantPairs <- localPairs
      importantTargets <- localTargets
      importantSources <- localSources
    }

    # Print results.
    if(verbose == TRUE){
      if(length(localPairs) < maxToPrint){
        print(paste("The pruned set of pairs is:", paste(sort(unlist(importantPairs)), collapse = ",")))
      }else{
        print(paste("The pruned set of pairs is:", paste(sort(unlist(importantPairs))[1:maxToPrint], collapse = ",")))
      }
    }
    return(list(pairs = unlist(importantPairs), targets = unlist(importantTargets), sources = unlist(importantSources)))
  })
  lengths <- unlist(lapply(prunedSubgraphs, function(g){return(length(g$pairs))}))
  prunedSubgraphs <- prunedSubgraphs[which(lengths > 0)]
  print(Sys.time() - tstart)
  return(prunedSubgraphs)
}

#' Obtain a prediction from a composite model, given the pairs to include in the model.
#' @param modelResults A ModelResults object.
#' @param trueVal The true values (predictions or outcomes) of the input data.
#' @param useActivation Whether or not to apply an activation function.
#' @return A vector of significance values named by model.
#' 
#' @export
ComputeIndividualPerformance <- function(modelResults, trueVal, pruningMethod, useActivation = FALSE){
  
  pred <- modelResults@model.input@edge.wise.prediction
  
  trueVal <- matrix(rep(trueVal, ncol(pred)), ncol = ncol(pred))
  tScore <- -1000

  # Compute the absolute error values.
  meanTrue <- matrix(rep(mean(trueVal), length(pred)), ncol = ncol(pred))
  meanAbsError <- abs(trueVal - meanTrue)
  predAbsError <- abs(trueVal - pred)
  
  # Compute Welch's t-test.
  n <- nrow(meanAbsError)
  meanMeanAbsError <- colMeans(meanAbsError)
  meanPredAbsError <- colMeans(predAbsError)
  meanMeanAbsErrorMat <- t(matrix(rep(meanMeanAbsError, n), ncol = n))
  meanPredAbsErrorMat <- t(matrix(rep(meanPredAbsError, n), ncol = n))
  stdevMeanAbsError <- sqrt(colSums((meanTrue - meanMeanAbsErrorMat) ^ 2) / n)
  stdevPredAbsError <- sqrt(colSums((pred - meanPredAbsErrorMat) ^ 2) / n)
  stdevMeanAbsErrorBar <- stdevMeanAbsError / sqrt(n)
  stdevPredAbsErrorBar <- stdevPredAbsError / sqrt(n)
  tScore <- (meanMeanAbsError - meanPredAbsError) / sqrt((stdevMeanAbsErrorBar ^ 2) + (stdevPredAbsErrorBar ^ 2))

  return(tScore)
}

#' Compute the weight of each predictor given the weights of different
#' metafeatures.
#' @param currentWeights The weights for each metafeature
#' @param metaFeatures A set of metafeatures, such as that found within ModelResults
#' @return A weight matrix for each sample and each predictor.
#' @export
ComputeMetaFeatureWeights <- function(metaFeatures, modelResults){
  weights <- lapply(1:length(metaFeatures), function(i){
    imp <- metaFeatures[[i]] * modelResults@current.metaFeature.weights[i]
    return(as.matrix(imp))
  })
  weights_all <- Reduce("+",weights)
  return(weights_all)
}
