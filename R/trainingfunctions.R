#' Run backpropagation for a single layer. In backpropagation, a gradient is
#' is computed by taking partial derivatives for each of the model weights. Given
#' a learning rate, the weights are adjusted according to the gradient.
#' @param modelResults An object of the ModelResults class.
#' @param prunedModels The models that remain after pruning.
#' @param Y.pred The predicted phenotype value.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
Backpropagate <- function(modelResults, prunedModels, Y.pred, minCutoff, maxCutoff, useCutoff = FALSE){
  # Calculate gradient.
  gradient <- computeGradient(modelResults = modelResults, prunedModels = prunedModels,
                              minCutoff = minCutoff,
                              maxCutoff = maxCutoff,
                              useCutoff = useCutoff)

  # Set previous weights.
  prevWeights <- modelResults@current.metaFeature.weights
  
  # Update gradient.
  modelResults@current.gradient <- as.matrix(gradient)
  modelResults@iteration.tracking[modelResults@current.iteration+1,
                                  which(grepl("Gradient", 
                                              colnames(modelResults@iteration.tracking)))] <- gradient
  
  # Set current weights, previous weights, and gradient.
  # Reference for the optimization algorithms: https://arxiv.org/abs/1609.04747
  if(modelResults@optimization.type == "SGD"){
    
    # Stochastic Gradient Descent
    modelResults@current.metaFeature.weights <- prevWeights - 
      modelResults@learning.rate * modelResults@current.gradient
    
    # Set the new weights and gradient.
    
  }else if(modelResults@optimization.type == "momentum"){
    # Momentum
    gamma <- 0.9
    update.vec <- gamma * modelResults@previous.update.vector + 
      modelResults@learning.rate * modelResults@current.gradient
    if(modelResults@current.iteration == 1){
      update.vec <- modelResults@learning.rate * modelResults@current.gradient
    }
    modelResults@current.metaFeature.weights <- prevWeights - update.vec
    modelResults@previous.update.vector <- update.vec
  }else if(modelResults@optimization.type == "adagrad"){
    # Adagrad
    G <- modelResults@sum.square.gradients
    epsilon <- 0.00000001
    update.vec <- (modelResults@learning.rate / sqrt(G+epsilon)) * modelResults@current.gradient
    if(modelResults@current.iteration == 1){
      update.vec <- modelResults@learning.rate * modelResults@current.gradient
    }
    modelResults@current.metaFeature.weights <- prevWeights - update.vec
    modelResults@sum.square.gradients <- modelResults@sum.square.gradients +
      (modelResults@current.gradient^2)
  }else if(modelResults@optimization.type == "adam"){
    # ADAM
    beta1 <- 0.9
    beta2 <- 0.999
    epsilon <- 0.00000001
    m <- beta1 * modelResults@previous.momentum + 
      (1-beta1) * modelResults@current.gradient
    v <- beta2 * modelResults@previous.update.vector +
      (1-beta2) * (modelResults@current.gradient^2)
    m.hat <- m / (1 - (beta1)^modelResults@current.iteration)
    v.hat <- v / (1 - (beta2)^modelResults@current.iteration)
    update.vec <- (modelResults@learning.rate * m.hat) / (sqrt(v.hat)+epsilon)
    modelResults@current.metaFeature.weights <- prevWeights - update.vec
    modelResults@previous.momentum <- m
    modelResults@previous.update.vector <- v
  }else if(modelResults@optimization.type == "newton"){
    # Calculate Hessian and update weights.
    hessian <- computeHessianSingleLayer(modelResults, convolution, pooling)
    update <- MASS::ginv(hessian) %*% modelResults@current.gradient
    modelResults@current.metaFeature.weights <- prevWeights - 
      (modelResults@learning.rate * update)
  }
  
  # Update the weights for this iteration.
  modelResults@iteration.tracking[modelResults@current.iteration+1,
                                  which(grepl("Weight", 
                                              colnames(modelResults@iteration.tracking)))] <- modelResults@current.metaFeature.weights
  
  # Return.
  return(modelResults)
}

#' Compute the Hessian for a single layer neural network.
#' Note that this function is currently not implemented for categorical outcomes.
#' @param modelResults An object of the ModelResults class.
#' @param convolution Whether or not to perform convolution
#' @param pooling Whether or not to perform pooling
computeHessianSingleLayer <- function(modelResults, convolution, pooling){
  # Components for derivative.
  A.hat <- modelResults@model.input@A.hat
  X <- modelResults@model.input@node.wise.prediction
  Theta.old <- matrix(rep(modelResults@current.metaFeature.weights, dim(X)[2]), ncol = dim(X)[2])
  Y <- modelResults@model.input@true.phenotypes
  
  # Convolve X.
  conv <- lapply(1:dim(X)[2], function(i){
    ret_val <- X[,i]
    if(convolution == TRUE){
      ret_val <- A.hat %*% X[,i]
    }
    return(ret_val)
  })
  
  # Compute filter.
  S <- modelResults@pooling.filter@filter
  S_all <- modelResults@pooling.filter@individual.filters
  
  # Compute the convolution values.
  convolution.val <- lapply(1:length(conv), function(i){
    return(conv[[i]])
  })
  if(pooling == TRUE){
    convolution.val <- lapply(1:length(S_all), function(i){
      S.flat <- rowSums(S_all[[i]])
      ret_val <- conv[[i]] * S.flat
      if(modelResults@weights.after.pooling == TRUE){
        ret_val <- t(conv[[i]]) %*% S_all[[i]]
      }
      return(ret_val)
    })
  }
  
  # Create array, where each matrix is one sample. These will be summed
  # together to obtain the Hessian.
  hessian.single.samples <- lapply(1:dim(X)[2], function(i){
    mat.X <- matrix(rep(X[,i], dim(X)[1]), nrow = dim(X)[1], ncol = dim(X)[1])
    return(mat.X * t(mat.X))
  })
  
  # Sum these together to obtain the Hessian.
  hessian <- hessian.single.samples[[1]]
  if(length(hessian.single.samples)>1){
    for(i in 2:length(hessian.single.samples)){
      hessian <- hessian + hessian.single.samples[[i]]
    }
  }
  
  # Return Hessian.
  return(hessian)
}

#' Compute the gradient for a single layer neural network.
#' @param modelResults An object of the ModelResults class.
#' @param prunedModels The models that remain after pruning.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
computeGradient <- function(modelResults, prunedModels, minCutoff, maxCutoff, useCutoff = FALSE){
  
  # Extract model components.
  S <- unlist(prunedModels)
  S_start <- unlist(lapply(S, function(pair){
    return(strsplit(pair, "__")[[1]][1])
  }))
  S_end <- unlist(lapply(S, function(pair){
    return(strsplit(pair, "__")[[1]][2])
  }))

  # Extract analyte types.
  analyteTypeOut <- modelResults@model.input@input.data@analyteType2
  analyteTypeIn <- modelResults@model.input@input.data@analyteType1
  if(S_start[1] %in% rownames(modelResults@model.input@input.data@analyteType2)){
    analyteTypeIn <- modelResults@model.input@input.data@analyteType2
  }
  if(S_end[1] %in% rownames(modelResults@model.input@input.data@analyteType1)){
    analyteTypeOut <- modelResults@model.input@input.data@analyteType1
  }
  Aind <- data.frame(analyteTypeIn[S_start,])
  Aout <- data.frame(analyteTypeOut[S_end,])

  # Extract covariates
  C <- do.call(cbind, lapply(modelResults@model.input@covariates, function(c){
    return(modelResults@model.input@input.data@sampleMetaData[,c])
  }))

  # Extract other components for gradient.
  P_hat <- CompositePrediction(pairs = S, modelResults = modelResults,
                               minCutoff = minCutoff,
                               maxCutoff = maxCutoff,
                               useCutoff = useCutoff)
  P <- modelResults@model.input@true.phenotypes
  Beta <- modelResults@model.input@model.properties[S,]
  M <- as.data.frame(do.call(rbind, modelResults@model.input@metaFeatures)[,S])
  rownames(M) <- names(modelResults@model.input@metaFeatures)
  phi <- modelResults@current.metaFeature.weights

  # Compute the gradient for each sample.
  sampGradients <- lapply(1:length(P), function(k){
    
    # Compute the derivative of P-hat.
    PhatDeriv <- 2 * (P[k] - P_hat[k])

    # Extract the independent, outcome, and covariate values.
    AindLocal <- data.frame(Aind[,k])
    if(length(AindLocal) == ncol(Aind) && ncol(Aind) > 1){
      AindLocal <- t(AindLocal)
      colnames(AindLocal) <- colnames(Aind)
    }
    AoutLocal <- data.frame(Aout[,k])
    if(length(AoutLocal) == ncol(Aout) && ncol(Aout) > 1){
      AoutLocal <- t(AoutLocal)
      colnames(AoutLocal) <- colnames(Aout)
    }
    CovarLocal <- NA
    if(length(modelResults@model.input@covariates) > 0){
      CovarLocal <- data.frame(C[k,])
      if(length(CovarLocal) == ncol(C)){
        CovarLocal <- t(CovarLocal)
        colnames(CovarLocal) <- colnames(C)
      }
    }

    # Compute the denominator.
    dhat <- Dhat(Aind = as.matrix(AindLocal), 
                 phi = phi, 
                 M = M, 
                 Beta2 = Beta[,"type"],
                 Beta3 = Beta[,"a:type"])
    
    # Compute the numerator.
    nhat <- Nhat(Aind = as.matrix(AindLocal),
                 Aout = as.matrix(AoutLocal),
                 C = CovarLocal,
                 phi = phi, 
                 M = M, 
                 Beta0 = Beta[,"(Intercept)"],
                 Beta1 = Beta[,"a"],
                 BetaC = Beta[,unlist(modelResults@model.input@covariates)])

    # Compute the derivative for each value of phi.
    phiGradients <- lapply(1:length(phi), function(gamma){
      # Extract the relevant component of M.
      Mgamma <- M[gamma,]

      # Compute the derivative of P-hat with respect to the denominator.
      dhatprime <- DhatPrime(Aind = as.matrix(AindLocal), 
                             M = Mgamma, 
                             Beta2 = Beta[,"type"],
                             Beta3 = Beta[,"a:type"])
      
      # Compute the derivative of P-hat with respect to the numerator.
      nhatprime <- NhatPrime(Aind = as.matrix(AindLocal),
                             Aout = as.matrix(AoutLocal),
                             C = CovarLocal,
                             M = Mgamma, 
                             Beta0 = Beta[,"(Intercept)"],
                             Beta1 = Beta[,"a"],
                             BetaC = Beta[,unlist(modelResults@model.input@covariates)])
      return((dhat * nhatprime - nhat * dhatprime) / (dhat ^ 2))
    })

    # Obtain the subgraph derivative.
    subgraphDerivDF <- as.data.frame(phiGradients)
    colnames(subgraphDerivDF) <- names(modelResults@model.input@metaFeatures)

    # If data are categorical, compute the derivative of the activation function.
    derivPredByAct <- 1
    if(modelResults@model.input@stype.class == "factor"){
      derivPredByAct <- DerivativeOfActivation(modelResults@activation.type, P_hat)
    }
    return(PhatDeriv * derivPredByAct * subgraphDerivDF)
  })
  sampGradientsDF <- do.call(rbind, sampGradients)
  gradients <- colSums(sampGradientsDF)
  
  # Return the derivative.
  return(-1 * gradients)
}

#' Compute the D-hat value, which is the denominator value of the composite prediction
#' for a sample k and a subgraph lambda.
#' @param Aind The matrix containing the values of the independent analyte type.
#' @param phi The weights of the meta-features.
#' @param M The meta-features.
#' @param Beta2 The phenotype term coefficients.
#' @param Beta3 The interaction term coefficients.
Dhat <- function(Aind, phi, M, Beta2, Beta3){
  phiMat <- t(matrix(rep(phi, nrow(Aind)), nrow = nrow(Aind)))
  sumWeights <- colSums(phiMat * M)
  term1 <- sum(Beta2 * sumWeights)
  term2 <- sum(Aind * Beta3 * sumWeights)
  return(term1 + term2)
}

#' Compute the D-hat value, which is the derivative of the denominator value of the composite prediction
#' for a sample k, a subgraph lambda, and a meta-feature metric gamma.
#' @param Aind The matrix containing the values of the independent analyte type.
#' @param M The meta-features.
#' @param Beta2 The phenotype term coefficients.
#' @param Beta3 The interaction term coefficients.
DhatPrime <- function(Aind, M, Beta2, Beta3){
  term1 <- sum(Beta2 * M)
  term2 <- sum(Aind * Beta3 * M)
  return(term1 + term2)
}

#' Compute the N-hat value, which is the numerator value of the composite prediction
#' for a sample k and a subgraph lambda.
#' @param Aind The matrix containing the values of the independent analyte type.
#' @param Aout The matrix containing the values of the outcome analyte type.
#' @param phi The weights of the meta-features.
#' @param M The meta-features.
#' @param Beta0 The intercept coefficients.
#' @param Beta1 The analyte coefficients.
#' @param BetaC The covariate coefficients.
#' @param C The covariates for the predictors of interest. If no covariates,
#' this value should be NA.
Nhat <- function(Aind, Aout, phi, M, Beta0, Beta1, BetaC, C){
  phiMat <- t(matrix(rep(phi, nrow(Aind)), nrow = nrow(Aind)))
  sumWeights <- colSums(phiMat * M)
  term1 <- sum(Aout * sumWeights)
  term2 <- sum(Beta0 * sumWeights)
  term3 <- sum(Aind * Beta1 * sumWeights)
  
  # Find term 4 if applicable.
  term4 <- 0
  if(!is.na(C)){
    covarMat <- matrix(rep(unlist(C), nrow(BetaC)), nrow = nrow(BetaC))
    sumWeightsMat <- matrix(rep(unlist(sumWeights), ncol(covarMat)), ncol = ncol(covarMat))
    term4 <- sum(BetaC * covarMat * sumWeightsMat)
  }
  return(term1 - (term2 + term3 + term4))
}

#' Compute the N-hat value, which is the derivative of the numerator value of the composite prediction
#' for a sample k, a subgraph lambda, and a meta-feature metric gamma.
#' @param Aind The matrix containing the values of the independent analyte type.
#' @param Aout The matrix containing the values of the outcome analyte type.
#' @param M The meta-feature values.
#' @param Beta0 The intercept coefficients.
#' @param Beta1 The analyte coefficients.
#' @param BetaC The covariate coefficients.
#' @param C The covariates for the predictors of interest. If no covariates,
#' this value should be NA.
NhatPrime <- function(Aind, Aout, M, Beta0, Beta1, BetaC, C){

  term1 <- sum(Aout * M)
  term2 <- sum(Beta0 * M)
  term3 <- sum(Aind * Beta1 * M)
  
  # Find term 4 if applicable.
  term4 <- 0
  if(!is.na(C)){
    covarMat <- matrix(rep(unlist(C), nrow(BetaC)), nrow = nrow(BetaC))
    sumWeightsMat <- matrix(rep(unlist(M), ncol(covarMat)), ncol = ncol(covarMat))
    term4 <- sum(BetaC * covarMat * sumWeightsMat)
  }
  
  return(term1 - (term2 + term3 + term4))
}

#' Compute the derivative of the error by the predictor.
#' @param activationType The name of the activation function.
#' @param pred The predicted value of a single sample for a single parameter
#' of interest.
DerivativeOfActivation <- function(activationType, pred){
  # Compute activation function.
  activation <- as.matrix(pred)
  d.act.d.pred <- as.matrix(rep(1, length(pred)))
  if(modelResults@model.input@outcome.type == "categorical"){
    if(modelResults@activation.type == "softmax"){
      activation <- SoftmaxWithCorrection(pred)
      # Set maximum and minimum values to prevent infinite and infinitecimal values.
      exp.pred <- exp(pred)
      exp.pred[which(exp.pred > 10000000)] <- 100000000
      exp.pred[which(exp.pred < 0.00000001)] <- 0.000000001
      d.act.d.pred <- (exp.pred * (exp.pred - sum(exp.pred)))/(sum(exp.pred)^2)
    }else if(modelResults@activation.type == "tanh"){
      activation <- TanhWithCorrection(pred)
      d.act.d.pred <- 2 * (1 - (tanh(2*(pred-1.5)))^2)
    }else if(modelResults@activation.type == "sigmoid"){
      activation <- SigmoidWithCorrection(pred)
      exp.pred <- exp(1-pred)
      exp.pred[which(exp.pred > 10000000)] <- 100000000
      exp.pred[which(exp.pred < 0.00000001)] <- 0.000000001
      d.act.d.pred <- -1 * exp.pred/((1+exp.pred)^2)
    }else if(modelResults@activation.type != "none"){
      stop(paste("Invalid activation type", modelResults@activation.type))
    }
  }
  
  # Return
  return(d.act.d.pred)
}

#' Apply the Softmax function to a numeric vector input. Note that we expect inputs
#' to be between 1 and 2, not 0 and 1. We therefore correct accordingly by subtracting
#' 1 from the input and adding 1 to the output.
#' @param input A vectorized numeric input.
SoftmaxWithCorrection <- function(input){
  numerator <- exp(input - 1)
  numerator[which(numerator > 10000000)] <- 100000000
  numerator[which(numerator < 0.00000001)] <- 0.000000001
  denominator <- sum(numerator)
  output <- (numerator/denominator)+1
  return(output)
}

#' Apply the Tanh function to a numeric vector input. Note that we expect inputs
#' to be between 1 and 2, not -1 and 1. We therefore correct accordingly.
#' @param input A vectorized numeric input.
TanhWithCorrection <- function(input){
  return(1.5 + tanh(2 * (input - 1.5)) / 2)
}

#' Apply the Sigmoid function to a numeric vector input. Note that we expect inputs
#' to be between 1 and 2, not -1 and 1. We therefore correct accordingly.
#' @param input A vectorized numeric input.
SigmoidWithCorrection <- function(input){
  exp.input <- exp(1-input)
  exp.input[which(exp.input > 10000000)] <- 100000000
  exp.input[which(exp.input < 0.00000001)] <- 0.000000001
  return(1 + (1 / (1 + exp.input)))
}

#' Predict Y given current weights.
#' @param modelResults An object of the ModelResults class.
#' @param prunedModels The models that remain after pruning. There should only
#' be one model at the end.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
DoPrediction <- function(modelResults, prunedModels, minCutoff = minCutoff, maxCutoff = maxCutoff, useCutoff = FALSE){
  
  # Obtain the final prediction.
  predictions <- CompositePrediction(pairs = unlist(prunedModels), modelResults = modelResults,
                                     minCutoff = minCutoff, 
                                     maxCutoff = maxCutoff, 
                                     useCutoff = FALSE)
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

#' Calculate minimum and maximum cutoffs for prediction.
#' @param modelResults An object of the ModelResults class.
CalculatePredictionCutoffs <- function(modelResults){
  fudgeFactor <- 0.1
  maxTrue <- max(modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype])
  minTrue <- min(modelResults@model.input@input.data@sampleMetaData[,modelResults@model.input@stype])
  marginMax <- fudgeFactor * abs(maxTrue)
  marginMin <- fudgeFactor * abs(minTrue)
  minCutoff <- minTrue - marginMin
  maxCutoff <- maxTrue + marginMax
  return(list(min = minCutoff, max = maxCutoff))
}

#' Train the graph learning model, using the specifications in the ModelResults
#' class and storing the results in the ModelResults class.
#' @param modelResults An object of the ModelResults class.
#' @param prunedModels The models that remain after pruning.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
DoSingleTrainingIteration <- function(modelResults, prunedModels, minCutoff, maxCutoff, useCutoff = FALSE){
  # Predict Y.
  Y.pred <- DoPrediction(modelResults, prunedModels,
                         minCutoff = minCutoff,
                         maxCutoff = maxCutoff,
                         useCutoff = useCutoff)
  modelResults@outcome.prediction <- as.numeric(as.matrix(Y.pred))

  # Backpropagate and calculate the error.
  modelResults <- Backpropagate(modelResults = modelResults,
                                prunedModels = prunedModels,
                                Y.pred = Y.pred,
                                minCutoff = minCutoff,
                                maxCutoff = maxCutoff,
                                useCutoff = useCutoff)

  # Modify the model results and return.
  return(modelResults)
}