#' Run backpropagation for a single layer. In backpropagation, a gradient is
#' is computed by taking partial derivatives for each of the model weights. Given
#' a learning rate, the weights are adjusted according to the gradient.
#' @param modelResults An object of the ModelResults class.
#' @param iteration The current iteration.
#' @param importance Importance metrics for each predictor and sample.
#' @param Y.pred The predicted phenotype value.
BackpropagateImportanceOnly <- function(modelResults, iteration, importance, Y.pred){
  # Calculate gradient.
  print(paste(Y.pred, modelResults@model.input@true.phenotypes))
  gradient <- computeGradientImportanceOnly(modelResults, importance, Y.pred)
  
  # Update gradient.
  modelResults@current.gradient <- as.matrix(gradient)
  modelResults@iteration.tracking[iteration+1,
                                  which(grepl("Gradient", 
                                              colnames(modelResults@iteration.tracking)))] <- gradient
  
  # Set current weights, previous weights, and gradient.
  # Reference for the optimization algorithms: https://arxiv.org/abs/1609.04747
  modelResults@previous.importance.weights <- modelResults@current.importance.weights
  if(modelResults@optimization.type == "SGD"){
    # Batch Gradient Descent
    modelResults@current.importance.weights <- modelResults@previous.importance.weights - 
      modelResults@learning.rate * modelResults@current.gradient
  }else if(modelResults@optimization.type == "momentum"){
    # Momentum
    gamma <- 0.9
    update.vec <- gamma * modelResults@previous.update.vector + 
      modelResults@learning.rate * modelResults@current.gradient
    if(modelResults@current.iteration == 1){
      update.vec <- modelResults@learning.rate * modelResults@current.gradient
    }
    modelResults@current.importance.weights <- modelResults@previous.importance.weights - update.vec
    modelResults@previous.update.vector <- update.vec
  }else if(modelResults@optimization.type == "adagrad"){
    # Adagrad
    G <- modelResults@sum.square.gradients
    epsilon <- 0.00000001
    update.vec <- (modelResults@learning.rate / sqrt(G+epsilon)) * modelResults@current.gradient
    if(modelResults@current.iteration == 1){
      update.vec <- modelResults@learning.rate * modelResults@current.gradient
    }
    modelResults@current.importance.weights <- modelResults@previous.importance.weights - update.vec
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
    modelResults@current.importance.weights <- modelResults@previous.importance.weights - update.vec
    modelResults@previous.momentum <- m
    modelResults@previous.update.vector <- v
  }else if(modelResults@optimization.type == "newton"){
    # Calculate Hessian and update weights.
    hessian <- computeHessianSingleLayer(modelResults, convolution, pooling)
    update <- MASS::ginv(hessian) %*% modelResults@current.gradient
    modelResults@current.importance.weights <- modelResults@previous.importance.weights - 
      (modelResults@learning.rate * update)
  }
  modelResults@iteration.tracking[iteration+1,
                                  which(grepl("Weight", 
                                              colnames(modelResults@iteration.tracking)))] <- modelResults@current.importance.weights
  
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
  Theta.old <- matrix(rep(modelResults@current.importance.weights, dim(X)[2]), ncol = dim(X)[2])
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
#' @param importance Importance metrics for each predictor and sample.
#' @param Y.pred The predicted phenotype value.
computeGradientImportanceOnly <- function(modelResults, importance, Y.pred){
  # Components for derivative.
  X <- modelResults@model.input@node.wise.prediction
  Y <- modelResults@model.input@true.phenotypes
  M <- importance
  
  # Compute the derivative for each sample.
  importanceWeightDeriv <- unlist(lapply(1:ncol(M),function(theta){
    
    # Compute the derivative for each sample. Remove the 2, which is a constant.
    derivErrByPred <- (Y - Y.pred) * Y.pred
    thetaSpecificPred <- sum(modelResults@current.importance.weights[theta] * 
                               M[theta] * Y.pred)
    derivPredByAct <- DerivativeOfActivation(modelResults@activation.type,
                                             thetaSpecificPred)
    derivActByTheta <- sum(M[theta] * Y.pred)
    
    # Return gradient.
    print(paste(derivErrByPred, derivPredByAct, derivActByTheta))
    return(derivErrByPred * derivPredByAct * derivActByTheta)
  }))

  # Return the derivative.
  return(importanceWeightDeriv)
}

#' Compute the derivative of the error by the predictor.
#' @param importance The importance values.
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