#' Sets up the test data so that it is in the correct format, then
#' runs prediction on the test data. To do this:
#' 1. Run pairwise prediction on the test data.
#' 2. Run pairwise prediction on the training data.
#' 3. Compute metafeatures on the test data.
#' 4. Predict the test values using the learned model.
#' @param inputDataTest An object of the IntLimData class corresponding to the test set.
#' @param model An object of the ModelResults class corresponding to the optimized model.
#' @param k The number of nearest neighbors to consider in localerr.
#' @param eigStep The number of eigenvectors to step by during the evaluation
#' in localerr.
#' Note that this must be less than the number of samples in localerr. Default = 10.
#' @param colIdInd The ID of the column that has the analyte ID's for the
#' independent variable. If blank, then the existing ID's are used.
#' @param colIdOut The ID of the column that has the analyte ID's for the
#' outcome variable. If blank, then the existing ID's are used.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @return A vector of final prediction values for the test data.
#' @export
DoTestSetupAndPrediction <- function(inputDataTest, model,
                                     k = 2, eigStep = 1,
                                     colIdInd = "databaseId",
                                     colIdOut = "databaseId",
                                     useCutoff = FALSE){
  # Get list of metafeatures.
  metaFeatureList <- names(model@model.input@metaFeatures)
  
  # Calculate prediction cutoffs.
  predictionCutoffs <- CalculatePredictionCutoffs(modelResults = model)
  
  # Run pairwise prediction for training data.
  pred.cv <- MultiOmicsGraphPrediction::RunPairwisePrediction(inputResults = model@model.input@model.properties,
                                                              inputData = model@model.input@input.data,
                                                              stype = model@model.input@stype,
                                                              covar = model@model.input@covariates,
                                                              independentVarType = model@model.input@independent.var.type,
                                                              outcomeType = model@model.input@outcome)
  pred.cv.test <- MultiOmicsGraphPrediction::RunPairwisePrediction(inputResults = model@model.input@model.properties,
                                                                   inputData = inputDataTest,
                                                                   stype = model@model.input@stype,
                                                                   covar = model@model.input@covariates,
                                                                   independentVarType = model@model.input@independent.var.type,
                                                                   outcomeType = model@model.input@outcome)
  
  # Get metafeatures for testing.
  metafeaturesTest <- MultiOmicsGraphPrediction::GetTestMetaFeatures(predictionsTrain = pred.cv,
                                                                     predictionsTest = pred.cv.test,
                                                                     stype = model@model.input@stype,
                                                                     inputDataTrain = model@model.input@input.data,
                                                                     inputDataTest = inputDataTest,
                                                                     metaFeatureList = metaFeatureList,
                                                                     k = k, eigStep = eigStep,
                                                                     modelStats = model@model.input@model.properties,
                                                                     colIdInd = colIdInd,
                                                                     colIdOut = colIdOut)
  
  # Predict testing values.
  weights <- ComputeMetaFeatureWeights(modelResults = model,
                                       metaFeatures = metafeaturesTest)
  predictions <- Predict(pairs = model@pairs,
                         inputData = inputDataTest, 
                         weights = weights, 
                         model = model,
                                useCutoff = useCutoff,
                                minCutoff = predictionCutoffs$min,
                                maxCutoff = predictionCutoffs$max,
                         outcomeType = model@model.input@outcome,
                         independentVarType = model@model.input@independent.var.type)
  return(predictions)
}

#' Splits input data into training and testing sets given a specified
#' set of testing samples. Alternatively, separate training and testing
#' sets can be read in prior to training. However, this function is a
#' useful helper for cross-validation.
#' @param inputData An object of the IntLimData class.
#' @param testingSamples Calculated importance metrics.
#' @return A list of two IntLimData objects: a training set and a testing set.
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
#' @param pairs A list of pairs to include in the composite model.
#' @param targets Target analytes for all pairs
#' @param sources Source analytes for all pairs
#' @param inputData The input testing data, which should be of an 
#' IntLimData class type.
#' @param weights The weights of each model, computed using ComputeMetaFeatureWeights()
#' @param model The learned model, an object of the modelResults class.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @param useActivation Whether or not to use the activation function if the phenotype
#' to predict is a factor. Default is TRUE. We set to FALSE during training.
#' @param independentVarType The independent variable type (1 or 2)
#' @param outcomeType The outcome type (1 or 2)
#' @return A vector of predictions
#' @export
Predict <- function(pairs, targets, sources, inputData, weights, model, minCutoff, maxCutoff, useCutoff = FALSE,
                    useActivation = TRUE, independentVarType, outcomeType){
  
  # Set variables for further analysis.
  covar <- model@model.input@covariates
  covariates <- model@model.input@model.properties
  analyteTgtVals <- inputData@analyteType1
  if(outcomeType == 2){
    analyteTgtVals <- inputData@analyteType2
  }
  analyteSrcVals <- inputData@analyteType2
  if(independentVarType == 1){
    analyteSrcVals <- inputData@analyteType1
  }
  covariateVals <- inputData@sampleMetaData
  wt <- weights[,pairs]
  covariatePairs <- covariates[pairs,]
  
  # The combination of these takes almost 2 seconds - why?
  # Analyte 1
  tgt <- t(analyteTgtVals[targets,])
  if(ncol(tgt) == nrow(wt) && nrow(tgt) == ncol(wt)){
    tgt <- analyteTgtVals[targets,]
  }
  weighted_a1 <- rowSums(wt * tgt)

  # beta0
  weighted_sum_b0 <- rowSums(wt * t(matrix(rep(covariatePairs[,"(Intercept)"], nrow(weights)), ncol = nrow(weights))))

  # beta1
  weighted_sum_b1 <- rowSums(wt * t(matrix(rep(covariatePairs[,"a"], nrow(weights)), ncol = nrow(weights))))

  # beta2
  weighted_sum_b2 <- rowSums(wt * t(matrix(rep(covariatePairs[,"type"], nrow(weights)), ncol = nrow(weights))))

  # beta3
  interactionTerm <- t(matrix(rep(covariatePairs[,"a:type"], nrow(weights)), ncol = nrow(weights)))
  src <- analyteSrcVals[sources,]
  if(ncol(src) == nrow(interactionTerm) && nrow(src) == ncol(interactionTerm)){
    src <- t(analyteSrcVals[sources,])
  }
  weighted_sum_b3 <- rowSums(wt * interactionTerm * src)

  # covariates
  weighted_sum_covars <- rep(0, nrow(weights))
  if(length(covar) > 0){
    weighted_sum_each <- lapply(covar, function(c){
      weighted_sum <- rowSums(wt * t(matrix(rep(covariatePairs[,c], nrow(weights)), ncol = nrow(weights)))
                                 * matrix(rep(covariateVals[,c], length(pairs)), ncol = length(pairs)))
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
  if(model@model.input@stype.class == "categorical" && useActivation == TRUE){
    if(model@activation.type == "softmax"){
      Y.pred <- round(SoftmaxWithCorrection(Y.pred))
    }else if(model@activation.type == "tanh"){
      Y.pred <- round(TanhWithCorrection(Y.pred))
    }else{
      Y.pred <- round(SigmoidWithCorrection(Y.pred))
    }
  }
  
  # Implement cutoff if desired.
  if(useCutoff == TRUE){
    Y.pred[which(Y.pred < minCutoff)] <- minCutoff
    Y.pred[which(Y.pred > maxCutoff)] <- maxCutoff
  }

  return(Y.pred)
}