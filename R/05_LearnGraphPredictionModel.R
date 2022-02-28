#' 
#' Format the input for graph-based learning. This input consists of:
#' 1. The Laplacian of a line graph built from the co-regulation graphs, where 
#' each node corresponds to a pair of analytes.
#' 2. A prediction value for each node of the line graph, for each sample X.
#' 3. The true prediction values Y for each sample X.
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param predictionGraphs A list of igraph objects, each of which includes
#' predictions for each edge.
#' @param coregulationGraph An igraph object containing the coregulation graph.
#' @param stype.class The class of the outcome ("numeric" or "categorical")
#' @param stype The phenotype or outcome of interest
#' @param edgeTypeList List containing one or more of the following to include
#' in the line graph:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
#' @param verbose Whether to print the number of predictions replaced in each sample.
#' TRUE or FALSE. Default is FALSE.
#' @export
FormatInput <- function(predictionGraphs, coregulationGraph,
                        inputData, stype.class, edgeTypeList, stype, verbose = FALSE){
  
  # Extract edge-wise predictions.
  predictions_by_node <- lapply(names(predictionGraphs), function(sampName){
    df_predictions <- igraph::as_data_frame(predictionGraphs[[sampName]])
    node_names <- paste(make.names(df_predictions$to), make.names(df_predictions$from),
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
  
  # Create a ModelInput object and return it.
  newModelInput <- methods::new("ModelInput", A.hat=A_hat, node.wise.prediction=t(predictions_flattened),
                                true.phenotypes=Y, outcome.type=stype.class, 
                                coregulation.graph=igraph::get.adjacency(coregulationGraph, sparse = FALSE), 
                                line.graph=as.matrix(A), modified.outliers=list())
  return(newModelInput)
}

#' Create the graph pooling filter, given the adjacency matrix of the input graph.
#' @param modelInputs An object of type "ModelInputs".
#' @param poolType One of "mean", "median", "max", or "min".
#' @export
CreatePoolingFilter <- function(modelInputs, poolType){
  
  # Convolve first.
  convolved <- modelInputs@A.hat %*% modelInputs@node.wise.prediction
  
  # Perform hierarchical clustering.
  hier <- IntLIM::doHierarchicalClustering(modelInputs)
  
  # Find the number of cuts to consider.
  heightRange <- c(min(hier$height), max(hier$height))
  heightRes <- 100
  cutSeq <- seq(from = heightRange[1], to = heightRange[2], 
                by = (heightRange[2] - heightRange[1])/heightRes)
  cutSeq <- cutSeq[order(-cutSeq)]
  
  # Find the maximum variance at each cut.
  bestHeight <- cutSeq[1]
  bestVar <- NULL
  for(c in 1:(length(cutSeq))){
    height <- cutSeq[c]
    cuts <- stats::cutree(hier, h = height)
    
    # For all cuts with more than one cluster, check the maximum variance
    # of the predictions.
    if(length(unique(cuts)) > 1){
      cutMaxVar <- lapply(unique(cuts), function(clust){
        currentClust <- names(cuts)[which(cuts == clust)]
        varReturn <- 0
        if(length(currentClust) > 1){
          cutVar <- unlist(lapply(1:ncol(convolved), function(j){
            return(var(convolved[currentClust,j]))
          }))
          varReturn <- max(cutVar)
        }
        return(varReturn)
      })
      mostlyBest <- quantile(unlist(cutMaxVar), 0.90)
      if(is.null(bestVar) || mostlyBest < bestVar){
        bestVar <- mostlyBest
        bestHeight <- c
      }
    }
  }
  
  # Arrange matrix of mappings.
  clusters <- stats::cutree(hier, h = height)[rownames(modelInputs@node.wise.prediction)]
  cluster_count <- length(unique(clusters))
  mappings <- matrix(0, ncol = cluster_count, nrow = nrow(convolved))
  cluster_length <- matrix(0, cluster_count)
  for(i in 1:cluster_count){
    which_in_cluster <- which(clusters == i)
    cluster_length[i] <- length(which_in_cluster)
    mappings[which_in_cluster, i] <- 1
  }
  newPoolingFilter <- methods::new("PoolingFilter", filter=mappings, filter.type=poolType,
                                   cluster.sizes=cluster_length, individual.filters=list())
  print(paste("Cutting the dendrogram at", bestHeight, "where the maximum predictor",
              "variance is", bestVar))
  return(newPoolingFilter)
}

#' Wrapper for CreatePoolingFilter.
#' @param modelInputs A list of objects of type "ModelInputs".
#' @param poolType One of "mean", "median", "max", or "min".
#' @export
CreatePoolingFilterAllFolds <- function(modelInputs, poolType){
  filter <- lapply(1:length(modelInputs), function(i){
    return(CreatePoolingFilter(modelInputs=modelInputs[[i]], poolType = poolType))
  })
  names(filter) <- paste("Fold", 1:length(modelInputs), sep = "_")
  return(filter)
}

#' Find edges that share nodes and add them to a data frame.
#' @param modelInputs An object of the ModelInput class.
#' @param poolingFilter A matrix that pools convolution results.
#' @param iterations Maximum number of iterations. Default is 1,000.
#' @param convergenceCutoff Cutoff for convergence. Default is 0.001.
#' @param stype.class One of either "character" or "numeric". Default is "numeric".
#' @param learningRate Learning rate to use during training. Default is 0.2
#' @param activationType Activation function. May be "softmax", "sigmoid", 
#' or "tanh". Default is "softmax". Irrelevant if stype.class is numeric.
#' @param gradientType One of either "stochastic" or "batch". Default is
#' "stochastic".
#' @param weightsAfterPooling Whether to include the weights after the pooling.
#' Default is FALSE.
#' @param optimizationType Type of optimization. May be "gradient", "momentum",
#' "adagrad", or "adam". Default is "adam".
#' @export
InitializeGraphLearningModel <- function(modelInputs, poolingFilter, 
                                         iterations = 1000,
                                         convergenceCutoff = 0.001, 
                                         stype.class = "numeric", 
                                         learningRate = 0.2,
                                         activationType = "softmax", 
                                         weightsAfterPooling = FALSE,
                                         gradientType = "stochastic",
                                         optimizationType = "adam"){
  
  # Initialize data frame with maximum number of iterations.
  weights_count <- dim(modelInputs@A.hat)[1]
  wt_name <- rownames(modelInputs@node.wise.prediction)
  if(weightsAfterPooling == TRUE){
    weights_count <- dim(poolingFilter@filter)[2]
    wt_name <- rep(1:dim(poolingFilter@filter)[2])
  }
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
  num_nodes <- nrow(modelInputs@node.wise.prediction)
  weights <- as.matrix(rep(1 / weights_count, weights_count))
  
  tracking.frame[1,3:(2+weights_count)] <- weights
  
  # Initialize and return results.
  newModelResults <- methods::new("ModelResults", model.input=modelInputs, pooling.filter=poolingFilter,
                        iteration.tracking=tracking.frame, max.iterations=iterations,
                        convergence.cutoff=convergenceCutoff, learning.rate=learningRate,
                        previous.weights=as.matrix(rep(0,length(weights))), 
                        current.weights=as.matrix(weights),
                        current.gradient=as.matrix(rep(-1,length(weights))),
                        previous.momentum=as.matrix(rep(0,length(weights))),
                        previous.update.vector=as.matrix(rep(0,length(weights))),
                        sum.square.gradients=as.matrix(rep(0,length(weights))),
                        current.iteration=0, activation.type=activationType,
                        weights.after.pooling=weightsAfterPooling,
                        optimization.type=optimizationType)
  return(newModelResults)
}

#' Wrapper for InitializeGraphLearningModel
#' @param modelInputs An object of the ModelInput class.
#' @param poolingFilter A matrix that pools convolution results.
#' @param iterations Maximum number of iterations. Default is 1,000.
#' @param convergenceCutoff Cutoff for convergence. Default is 0.001.
#' @param stype.class One of either "character" or "numeric". Default is "numeric".
#' @param learningRate Learning rate to use during training. Default is 0.2
#' @param activationType Activation function. May be "softmax", "sigmoid", 
#' or "tanh". Default is "softmax". Irrelevant if stype.class is numeric.
#' @param gradientType One of either "stochastic" or "batch". Default is
#' "stochastic".
#' @param weightsAfterPooling Whether to include the weights after the pooling.
#' Default is FALSE.
#' @param optimizationType Type of optimization. May be "gradient", "momentum",
#' "adagrad", or "adam". Default is "adam".
#' @export
InitializeGraphLearningModelAllFolds <- function(modelInputs, poolingFilter, 
                                              iterations = 1000,
                                              convergenceCutoff = 0.001, 
                                              stype.class = "numeric", 
                                              learningRate = 0.2,
                                              activationType = "softmax", 
                                              weightsAfterPooling = FALSE,
                                              gradientType = "stochastic",
                                              optimizationType = "adam"){
  filter <- lapply(1:length(modelInputs), function(i){
    return(InitializeGraphLearningModel(modelInputs = modelInputs[[i]], 
                                        poolingFilter = poolingFilter[[i]], 
                                        iterations = iterations,
                                        convergenceCutoff = convergenceCutoff, 
                                        stype.class = stype.class, 
                                        learningRate = learningRate,
                                        activationType = activationType, 
                                        weightsAfterPooling = weightsAfterPooling,
                                        gradientType = gradientType,
                                        optimizationType = optimizationType))
  })
  names(filter) <- paste("Fold", 1:length(modelInputs), sep = "_")
  return(filter)
}

#' Run least squares optimization to optimize the weights. This function
#' assumes that the weights are applied after pooling.
#' @param modelResults An object of the ModelResults class.
#' @export
RunLeastSquaresOptimization <- function(modelResults){
  # Get Components.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Y.pred <- X
  
  # Find Median Pools.
  S_all <- AdjustFilter(poolingFilter = modelResults@pooling.filter,
                        A.hat = A.hat, 
                        X = X)
  
  # Convolve and pool.
  Y.pred.list <- lapply(1:dim(X)[2], function(i){
    return(t(A.hat %*% X[,i])  %*% S_all[[i]])
  })
  Y.pred <- do.call(rbind, Y.pred.list)
  
  # Train.
  model <- glmnet::glmnet(x = Y.pred, y = modelResults@model.input@true.phenotypes,
                  intercept = FALSE, lower.limits = 0, lambda = 0)
  
  # Return.
  return(model)
}

#' Wrapper for RunLeastSquaresOptimization.
#' @param modelResults A list of objects of the ModelResults class.
#' @export
RunLeastSquaresOptimizationAllFolds <- function(modelResults){
  result <- lapply(1:length(modelResults), function(i){
    return(RunLeastSquaresOptimization(modelResults=modelResults[[i]]))
  })
  names(result) <- paste("Fold", 1:length(modelResults), sep = "_")
  return(result)
}

#' Predict test data using a least squares model. This function
#' assumes that the weights are applied after pooling.
#' @param modelInput An object of the ModelInput class. This input corresponds
#' to the test data.
#' @param trainingModelResults An object of the ModelResults class. This corresponds to
#' the modelResults object for the training data.
#' @param leastSquaresModel A least squares model learned using 
#' RunLeastSquaresOptimization.
#' @export
PredictFromLeastSquares <- function(modelInput, trainingModelResults, leastSquaresModel){
  # Construct testing input.
  Y.pred.test <- modelInput@node.wise.prediction
  A.hat = trainingModelResults@model.input@A.hat
  
  # Adjust filter.
  S_all <- AdjustFilter(poolingFilter = trainingModelResults@pooling.filter, 
                        A.hat = A.hat,
                        X = Y.pred.test)
  
  # Convolve and pool.
  Y.pred.list.test <- lapply(1:dim(Y.pred.test)[2], function(i){
    return(t(A.hat %*% Y.pred.test[,i])  %*% S_all[[i]])
  })
  Y.pred.test <- do.call(rbind, Y.pred.list.test)
  test.predict <- stats::predict(leastSquaresModel, Y.pred.test)
  rownames(test.predict) <- colnames(modelInput@node.wise.prediction)
  return(test.predict)
}

#' Wrapper for RunLeastSquaresOptimization.
#' @param modelInput A list of objects of the ModelInput class. This input corresponds
#' to the test data.
#' @param trainingModelResults A list of objects of the ModelResults class. This corresponds to
#' the modelResults object for the training data.
#' @param leastSquaresModel A list of least squares models learned using 
#' RunLeastSquaresOptimization.
#' @export
PredictFromLeastSquaresAllFolds <- function(modelInput, trainingModelResults, leastSquaresModel){
  result <- lapply(1:length(modelInput), function(i){
    return(PredictFromLeastSquares(modelInput = modelInput[[i]],
                                   trainingModelResults = trainingModelResults[[i]],
                                   leastSquaresModel = leastSquaresModel[[i]]))
  })
  names(result) <- paste("Fold", 1:length(modelInput), sep = "_")
  return(result)
}

#' Train the graph learning model, using the specifications in the ModelResults.
#' class and storing the results in the ModelResults class.
#' @param modelResults An object of the ModelResults class.
#' @param pooling Whether or not to perform pooling during training.
#' @param convolution Whether or not to perform convolution during training.
#' @param stochastic Whether or not training should be stochastic.
#' @param ridgeRegressionWeight The hyperparameter weight assigned
#' to the ridge regression parameter (often referred to as lambda in the
#' literature)
#' @param varianceWeight The hyperparameter weight assigned to the difference
#' in variances between Y and the predicted value of Y.
#' @param verbose Whether to print results as you run the model.
#' @export
TrainGraphLearningModel <- function(modelResults, pooling, convolution,
                                    stochastic, ridgeRegressionWeight,
                                    varianceWeight, verbose = TRUE){
  
  # Start the first iteration and calculate a dummy weight delta.
  modelResults@current.iteration <- 1
  modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
    modelResults@current.iteration
  weight.delta <- sqrt(sum((modelResults@current.weights - modelResults@previous.weights)^2))
  
  # Repeat the training process for all iterations, until the maximum is reached
  # or until convergence.
  while(modelResults@current.iteration < (modelResults@max.iterations - 1)
        && (weight.delta > modelResults@convergence.cutoff)){
    # For stochastic training, permute the samples, then compute the gradient one
    # sample at a time.
    # For batch training, compute the gradient over all samples.
    if(stochastic == TRUE){
      # Permute samples.
      perm_samples <- sample(1:dim(modelResults@model.input@node.wise.prediction)[2],
                             dim(modelResults@model.input@node.wise.prediction)[2])
      for(i in perm_samples){
        # Do training iteration for each sample.
        newModelResults <- modelResults
        newModelResults@model.input@node.wise.prediction <- 
          as.matrix(modelResults@model.input@node.wise.prediction[,i])
        newModelResults@model.input@true.phenotypes <- 
          modelResults@model.input@true.phenotypes[i]
        newModelResults <- DoSingleTrainingIteration(newModelResults, modelResults@current.iteration,
                                                  pooling, convolution, ridgeRegressionWeight,
                                                  varianceWeight)
        # Update weights and gradient in the model results according to the
        # results of this sample.
        modelResults@current.weights <- newModelResults@current.weights
        modelResults@previous.weights <- newModelResults@previous.weights
        modelResults@current.gradient <- newModelResults@current.gradient
        modelResults@outcome.prediction <- newModelResults@outcome.prediction
        modelResults@previous.momentum <- newModelResults@previous.momentum
        modelResults@previous.update.vector <- newModelResults@previous.update.vector
        modelResults@iteration.tracking <- newModelResults@iteration.tracking
      }
      # Compute the prediction error over all samples.
      Y.pred <- DoPrediction(modelResults, modelResults@current.iteration,
                                                           pooling, convolution)
      if(modelResults@model.input@outcome.type == "categorical"){
        modelResults@iteration.tracking$Error[modelResults@current.iteration+1] <- 
          ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
      }else{
        modelResults@iteration.tracking$Error[modelResults@current.iteration+1] <- 
          ComputeNRMSE(modelResults@model.input@true.phenotypes, Y.pred)
      }
      
    # This is the batch case.  
    }else{
      modelResults <- DoSingleTrainingIteration(modelResults, modelResults@current.iteration,
                                                pooling, convolution, ridgeRegressionWeight,
                                                varianceWeight)
    }
    
    # Update the iteration.
    modelResults@current.iteration <- modelResults@current.iteration + 1
    modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
      modelResults@current.iteration
    
    # Print the weight delta and error.
    weight.delta <- sqrt(sum((modelResults@current.weights - modelResults@previous.weights)^2))
    if(modelResults@current.iteration %% 1 == 0 && verbose == TRUE){
      print(paste("iteration", modelResults@current.iteration, ": weight delta is", weight.delta,
                  "and error is", 
                  modelResults@iteration.tracking$Error[modelResults@current.iteration]))
    }
  }
  # If we exited before the maximum number of iterations, remove the rest of the
  # tracking data.
  if(modelResults@current.iteration < modelResults@max.iterations){
    modelResults@iteration.tracking <- modelResults@iteration.tracking[1:modelResults@current.iteration,]
  }
  return(modelResults)
}

#' Run a prediction on new data using the graph learning model, and compute
#' the absolute error delta values.
#' @param weights The list of all learned weights.
#' @param convolutions Number of convolutions to perform. This should be the
#' same number of convolutions used in learning the weights.
#' @param modelInput A list of objects of the ModelInput class.
#' @param minimum The minimum prediction value to allow. Default is NULL.
#' @param maximum The maximum prediction value to allow. Default is NULL.
#' @export
GetErrorDeltas <- function(weights, modelInput, convolutions, minimum = NULL, 
                           maximum = NULL){
  # Convolve the input.
  Y.formula <- modelInput@node.wise.prediction
  if(convolutions > 0){
    
    # Modify the convolutional matrix if more than
    # one convolution is desired.
    if(convolutions > 1){
      for(c in 2:convolutions){
        modelInput@A.hat <- modelInput@A.hat %*% modelInput@A.hat
      }
    }
    
    # Do convolution.
    Y.formula <- modelInput@A.hat %*% Y.formula
  }
  deltas <- unlist(lapply(1:length(modelInput@true.phenotypes),function(j){
    solution <- sum(Y.formula[,j] * weights)
    # Adjust to fit minimum and maximum.
    if(!is.null(minimum)){
      solution[which(solution < minimum)] <- minimum
    }
    if(!is.null(maximum)){
      solution[which(solution > maximum)] <- maximum
    }  
    phenotype <- modelInput@true.phenotypes[j]
    return(unlist(unname(abs(solution - phenotype))))
  }))
  names(deltas) <- names(modelInput@true.phenotypes)
  return(deltas)
}

#' Run a prediction on new data using the graph learning model, and compute
#' the absolute error delta values.
#' @param weights The list of all learned weights.
#' @param convolutions Number of convolutions to perform. This should be the
#' same number of convolutions used in learning the weights.
#' @param modelInput A list of objects of the ModelInput class.
#' @param minimum The minimum prediction value to allow. Default is NULL.
#' @param maximum The maximum prediction value to allow. Default is NULL.
#' @export
GetErrorDeltasAllFolds <- function(weights, modelInput, convolutions, minimum = NULL, 
                                   maximum = NULL){
  return(lapply(1:length(modelInput), function(i){
    # Get error deltas.
    return(GetErrorDeltas(weights=weights[[i]], modelInput=modelInput[[i]], 
                          convolutions=convolutions, minimum=minimum,
                          maximum=maximum))
  }))
}
#' Run a prediction on new data using the graph learning model.
#' @param modelResults An object of the ModelResults class.
#' @param pooling Whether or not to pool the weights.
#' @param convolution Whether or not to perform convolution.
#' @param testInput An object of the ModelInput class.
#' @export
PredictTesting <- function(modelResults, pooling, convolution, testInput){
  # Propagate forward.
  A.hat <- modelResults@model.input@A.hat
  X <- testInput@node.wise.prediction
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
  
  # Calculate error.
  error <- 1
  if(modelResults@model.input@outcome.type == "categorical"){
    error <- ComputeClassificationError(testInput@true.phenotypes, Y.pred)
  }else{
    error <- ComputeNRMSE(testInput@true.phenotypes, Y.pred)
  }
  
  # Modify the model results and return.
  return(list("Y.pred" = Y.pred, "Error" = error))
}

#' Compute classification error.
#' @param true.Y The true phenotype of each sample.
#' @param pred.Y The predicted phenotype of each sample.
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
#' @export
ComputeNRMSE <- function(true.Y, pred.Y){
  RMSD <- sqrt(sum((true.Y - pred.Y)^2) / length(true.Y))
  NRMSE <- RMSD
  if(length(true.Y) > 1){
    NRMSE <- RMSD / (max(true.Y) - min(true.Y))
  }
  return(NRMSE)
}