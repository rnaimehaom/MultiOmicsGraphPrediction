#' Find edges that share nodes and add them to a data frame.
#' @param modelInputs An object of the ModelInput class.
#' @param importance Calculated importance metrics.
#' @param iterations Maximum number of iterations. Default is 1,000.
#' @param convergenceCutoff Cutoff for convergence. Default is 0.001.
#' @param modelType The type of model the user wishes to learn. Users may
#' select from the following:
#' - "importance": A linear combination of predictors, where each predictor's
#'   weight is learned as a function of importance metrics.
#' @param stype.class One of either "character" or "numeric". Default is "numeric".
#' @param learningRate Learning rate to use during training. Default is 0.2
#' @param activationType Activation function. May be "softmax", "sigmoid", 
#' "tanh", or "none". Default is "none", meaning stype.class is numeric.
#' @param optimizationType Type of optimization. May be "SGD", "momentum",
#' "adagrad", or "adam". Default is "SGD".
#' @param initialImportanceWeights Initial weights for model importance. Default
#' is NULL, which results in each importance metric being given equal weight.
#' @export
InitializeGraphLearningModel <- function(modelInputs,
                                         iterations = 1000,
                                         convergenceCutoff = 0.0000000001,
                                         modelType = "importance",
                                         stype.class = "numeric", 
                                         learningRate = 0.2,
                                         activationType = "none", 
                                         optimizationType = "SGD",
                                         initialImportanceWeights = NULL){
  # Initialize importance weights.
  weights_count <- length(modelInputs@importance)
  wt_name <- names(modelInputs@importance)
  
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
  num_nodes <- nrow(modelInputs@node.wise.prediction)
  weights <- as.matrix(rep(1 / weights_count, weights_count))
  if(!is.null(initialImportanceWeights)){
    weights <- initialImportanceWeights
  }
  
  tracking.frame[1,3:(2+weights_count)] <- weights

  # Initialize and return results.
  newModelResults <- methods::new("ModelResults", model.input=modelInputs,
                        iteration.tracking=tracking.frame, max.iterations=iterations,
                        convergence.cutoff=convergenceCutoff, learning.rate=learningRate,
                        previous.importance.weights=as.matrix(rep(0,length(weights))), 
                        current.importance.weights=as.matrix(weights),
                        current.gradient=as.matrix(rep(-1,length(weights))),
                        previous.momentum=as.matrix(rep(0,length(weights))),
                        previous.update.vector=as.matrix(rep(0,length(weights))),
                        sum.square.gradients=as.matrix(rep(0,length(weights))),
                        current.iteration=0, activation.type=activationType,
                        optimization.type=optimizationType)
  return(newModelResults)
}

#' Optimize the combination of predictors by importance alone (in other words,
#' exclude pooling and combine in a single layer using a linear combination).
#' @param modelResults An object of the ModelResults class.
#' @param verbose Whether to print results as you run the model.
#' @export
OptimizeImportanceCombo <- function(modelResults, verbose = TRUE){
  
  # Start the first iteration and calculate a dummy weight delta.
  modelResults@current.iteration <- 1
  modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
    modelResults@current.iteration
  weight.delta <- sqrt(sum((modelResults@current.importance.weights - 
                              modelResults@previous.importance.weights)^2))
  
  # Get the initial pruned models.
  pairsPredAll <- MultiOmicsGraphPrediction::ObtainSubgraphNeighborhoods(modelInput = modelResults@model.input, percentOverlapCutoff = 50)
  prunedModels <- MultiOmicsGraphPrediction::DoSignificancePropagation(pairs = pairsPredAll, modelResults = modelResults,
                                                                          verbose = FALSE, makePlots = FALSE)

  # Placeholder for predictions.
  Y.pred <- rep(0,nrow(modelResults@model.input@node.wise.prediction))
  
  # Repeat the training process for all iterations, until the maximum is reached
  # or until convergence.
  str(modelResults@convergence.cutoff)
  while(modelResults@current.iteration < (modelResults@max.iterations - 1)
        && (weight.delta > modelResults@convergence.cutoff)){
    # For stochastic training, permute the samples, then compute the gradient one
    # sample at a time.
    perm_samples <- sample(1:nrow(modelResults@model.input@node.wise.prediction),
                           nrow(modelResults@model.input@node.wise.prediction))

    for(i in perm_samples){
      # Do training iteration for each sample.
      newModelResults <- modelResults
      newModelResults@model.input@node.wise.prediction <- 
        as.matrix(modelResults@model.input@node.wise.prediction[i,])
      newModelResults@model.input@true.phenotypes <- 
        modelResults@model.input@true.phenotypes[i]
      importanceSamp <- lapply(1:length(modelResults@model.input@importance), function(imp){
        df <- t(as.data.frame(modelResults@model.input@importance[[imp]][i,]))
        rownames(df) <- rownames(modelResults@model.input@node.wise.prediction)[i]
        colnames(df) <- colnames(modelResults@model.input@node.wise.prediction)
        return(df)
      })
      names(importanceSamp) <- names(modelResults@model.input@importance)
      newModelResults@model.input@importance <- importanceSamp
      newModelResults <- DoSingleTrainingIteration(modelResults = newModelResults,
                                                   prunedModels = prunedModels)
      
      # Update weights and gradient in the model results according to the
      # results of this sample.
      modelResults@current.importance.weights <- newModelResults@current.importance.weights
      modelResults@previous.importance.weights <- newModelResults@previous.importance.weights
      modelResults@current.gradient <- newModelResults@current.gradient
      modelResults@outcome.prediction <- newModelResults@outcome.prediction
      modelResults@previous.momentum <- newModelResults@previous.momentum
      modelResults@previous.update.vector <- newModelResults@previous.update.vector
      modelResults@iteration.tracking <- newModelResults@iteration.tracking
      
      # Fill in the prediction for error calculation.
      Y.pred <- DoPrediction(modelResults = newModelResults, prunedModels = prunedModels)
    }
    # Compute the prediction error over all samples.
    if(modelResults@model.input@outcome.type == "categorical"){
      modelResults@iteration.tracking$Error[modelResults@current.iteration] <- 
        ComputeClassificationError(modelResults@model.input@true.phenotypes, Y.pred)
    }else{
      modelResults@iteration.tracking$Error[modelResults@current.iteration] <- 
        ComputeRMSE(modelResults@model.input@true.phenotypes, Y.pred)
    }
    currentError <- modelResults@iteration.tracking$Error[modelResults@current.iteration]
    
    # Update the iteration.
    modelResults@current.iteration <- modelResults@current.iteration + 1
    modelResults@iteration.tracking$Iteration[modelResults@current.iteration+1]<-
      modelResults@current.iteration
    
    # Get the new pruned models.
    prunedModels <- MultiOmicsGraphPrediction::DoSignificancePropagation(pairs = pairsPredAll, modelResults = modelResults,
                                                                         verbose = FALSE, makePlots = FALSE)
    prunedModelSizes <- lapply(prunedModels, function(model){return(length(model))})
    
    # Print the weight delta and error.
    weight.delta <- sqrt(sum((modelResults@current.importance.weights - modelResults@previous.importance.weights)^2))
    if(modelResults@current.iteration %% 1 == 0 && verbose == TRUE){
      print(paste("iteration", modelResults@current.iteration, ": weight delta is", weight.delta,
                  "and error is", paste0(currentError, ". Subgraph sizes are ", 
                  paste(prunedModelSizes, collapse = ", "))))
    }
  }
  # If we exited before the maximum number of iterations, remove the rest of the
  # tracking data.
  if(modelResults@current.iteration < modelResults@max.iterations){
    modelResults@iteration.tracking <- modelResults@iteration.tracking[1:modelResults@current.iteration,]
  }
  return(modelResults)
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
  S_all <- lapply(1:dim(X)[2], function(i){
    return(S)
  })
  
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
  
  # Adjust filter. Need to fill this in.
  S_all <- lapply(1:dim(X)[2], function(i){
    return(NULL)
  })
  
  # Convolve and pool.
  Y.pred.list.test <- lapply(1:dim(Y.pred.test)[2], function(i){
    return(t(A.hat %*% Y.pred.test[,i])  %*% S_all[[i]])
  })
  Y.pred.test <- do.call(rbind, Y.pred.list.test)
  test.predict <- stats::predict(leastSquaresModel, Y.pred.test)
  rownames(test.predict) <- colnames(modelInput@node.wise.prediction)
  return(test.predict)
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
  Theta.old <- matrix(rep(modelResults@current.importance.weights, dim(X)[2]), ncol = dim(X)[2])
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
        # Need to fill this in.
        S_all <- NULL#CreateFilter(poolingFilter = modelResults@pooling.filter, A.hat = A.hat, 
                              #X = X)
        modelResults@pooling.filter@individual.filters <- S_all
      }
      Y.pred <- unlist(lapply(1:length(S_all), function(i){
        return(sum(t(Y.pred[,i]) %*% S_all[[i]] * Theta.old[,i]))
      }))
    }else{
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
    error <- ComputeRMSE(testInput@true.phenotypes, Y.pred)
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
ComputeRMSE <- function(true.Y, pred.Y){
  RMSD <- sqrt(sum((true.Y - pred.Y)^2) / length(true.Y))
  return(RMSD)
}