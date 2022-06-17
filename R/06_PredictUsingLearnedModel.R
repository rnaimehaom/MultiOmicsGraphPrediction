#' Splits input data into training and testing sets given a specified
#' set of testing samples. Alternatively, separate training and testing
#' sets can be read in prior to training. However, this function is a
#' useful helper for cross-validation.
#' @param inputData An object of the IntLimData class.
#' @param testingSamples Calculated importance metrics.
#' @export
SplitTrainingAndTesting <- function(inputData, testingSamples){
  
  # Initialize training and testing data to be identical to original data sets.
  inputDataTrain <- inputData
  inputDataTest <- inputData
  
  # Remove testing sample from training data.
  trainingSamples <- setdiff(rownames(inputDataTrain@sampleMetaData), testingSamples)
  inputDataTrain@analyteType1 <- as.matrix(inputDataTrain@analyteType1[,trainingSamples])
  colnames(inputDataTrain@analyteType1) <- trainingSamples
  inputDataTrain@analyteType2 <- inputDataTrain@analyteType2[,trainingSamples]
  colnames(inputDataTrain@analyteType2) <- trainingSamples
  inputDataTrain@sampleMetaData <- inputDataTrain@sampleMetaData[trainingSamples,]
  rownames(inputDataTrain@sampleMetaData) <- trainingSamples
  
  # Retain only testing sample for testing data.
  inputDataTest@analyteType1 <- as.matrix(inputDataTest@analyteType1[,testingSamples])
  colnames(inputDataTest@analyteType1) <- testingSamples
  inputDataTest@analyteType2 <- as.matrix(inputDataTest@analyteType2[,testingSamples])
  colnames(inputDataTest@analyteType2) <- testingSamples
  inputDataTest@sampleMetaData <- as.data.frame(inputDataTest@sampleMetaData[testingSamples,])
  rownames(inputDataTest@sampleMetaData) <- testingSamples
  
  return(list(training = inputDataTrain, testing = inputDataTest))
}

#' Run a prediction on new data using the graph learning model.
#' @param inputData The input testing data, which should be of an 
#' IntLimData class type.
#' @param metafeatures The metafeatures calculated using GetTestMetaFeatures()
#' @param model The learned model, an object of the modelResults class.
#' @export
PredictTesting <- function(inputData, metafeatures, model){
  
  # Set variables for further analysis.
  covar <- model@model.input@covariates
  covariates <- model@model.input@model.properties
  analyteTgtVals <- inputData@analyteType1
  if(model@model.input@outcome == 2){
    analyteTgtVals <- inputData@analyteType2
  }
  analyteSrcVals <- inputData@analyteType2
  if(model@model.input@independent.var.type == 1){
    analyteSrcVals <- inputData@analyteType1
  }
  covariateVals <- inputData@sampleMetaData
  weights <- ComputeMetaFeatureWeights(modelResults = model,
                                       metaFeatures = metafeatures)
  
  # Analyte 1
  weighted_a1 <- rep(0, nrow(weights))
  for(i in 1:length(model@pairs)){
    pair <- model@pairs[[i]]
    weighted_a1 <- weighted_a1 + weights[,pair] * analyteTgtVals[strsplit(pair, "__")[[1]][2],]
  }
  
  # beta0
  weighted_sum_b0 <- rep(0, nrow(weights))
  for(i in 1:length(model@pairs)){
    pair <- model@pairs[[i]]
    weighted_sum_b0 <- weighted_sum_b0 + weights[,pair] * rep(covariates[pair,"(Intercept)"], nrow(weights))
  }
  
  # beta1
  weighted_sum_b1 <- rep(0, nrow(weights))
  for(i in 1:length(model@pairs)){
    pair <- model@pairs[[i]]
    weighted_sum_b1 <- weighted_sum_b1 + weights[,pair] * rep(covariates[pair, "a"], nrow(weights)) * 
      analyteSrcVals[strsplit(pair, "__")[[1]][1],]
  }
  
  # beta2
  weighted_sum_b2 <- rep(0, nrow(weights))
  for(i in 1:length(model@pairs)){
    pair <- model@pairs[[i]]
    weighted_sum_b2 <- weighted_sum_b2 + weights[,pair] * rep(covariates[pair, "type"], nrow(weights))
  }
  
  # beta3
  weighted_sum_b3 <- rep(0, nrow(weights))
  for(i in 1:length(model@pairs)){
    pair <- model@pairs[[i]]
    weighted_sum_b3 <- weighted_sum_b3 + weights[,pair] * rep(covariates[pair, "a:type"], nrow(weights))* 
      analyteSrcVals[strsplit(pair, "__")[[1]][1],]
  }
  
  # covariates
  weighted_sum_covars <- rep(0, nrow(weights))
  if(length(covar) > 0){
    weighted_sum_each <- lapply(covar, function(c){
      weighted_sum <- rep(0, nrow(weights))
      for(i in 1:length(model@pairs)){
        pair <- model@pairs[[i]]
        tryCatch({
          weighted_sum <- weighted_sum + weights[,pair] * rep(covariates[pair, c], nrow(weights)) *
            covariateVals[,c]
        }, error = function(e){
          print(e)
          print(weighted_sum)
          print(weights[,pair])
          print(rep(covariates[pair, c], nrow(weights)))
          print(covariateVals[,c])
        })
        
      }
      return(weighted_sum)
    })
    weighted_sum_covars <- Reduce('+', weighted_sum_each)
  }
  
  # Final value.
  
  denom <- weighted_sum_b2 + weighted_sum_b3
  denom[which(denom == 0)] <- 0.0001
  Y.pred <- (weighted_a1 - weighted_sum_b0 - weighted_sum_b1 - weighted_sum_covars) / denom
  
  
  # Use activation function if output is of a character type. Note that all character
  # types are converted into factors, and since only binary factors are accepted by
  # the package, the values will be 1 (for the alphanumerically lowest level) and 2
  # (for the alphanumerically highest level).
  if(model@model.input@outcome.type == "categorical"){
    if(model@activation.type == "softmax"){
      Y.pred <- round(SoftmaxWithCorrection(Y.pred))
    }else if(model@activation.type == "tanh"){
      Y.pred <- round(TanhWithCorrection(Y.pred))
    }else{
      Y.pred <- round(SigmoidWithCorrection(Y.pred))
    }
  }
  
  return(Y.pred)
}