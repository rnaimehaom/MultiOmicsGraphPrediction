#' Sets up the test data so that it is in the correct format, then
#' runs prediction on the test data. To do this:
#' 1. Run pairwise prediction on the test data.
#' 2. Run pairwise prediction on the training data.
#' 3. Compute metafeatures on the test data.
#' 4. Predict the test values using the learned model.
#' @param inputDataTrain An object of the IntLimData class corresponding to the training set.
#' @param inputDataTest An object of the IntLimData class corresponding to the test set.
#' @param model An object of the ModelResults class corresponding to the optimized model.
#' @param stype The variable to predict.
#' @param covar The clinical covariates to include in the model. These should be the same
#' covariates that were included when running the IntLIM linear models.
#' @param independentVarType The independent variable type (1 or 2)
#' @param outcomeType The outcome type (1 or 2)
#' @param metaFeatureList A list of the valid metrics to include. Valid metrics are
#' "pdf", "localerr", "globalerr", and "pathway".
#' @param k The number of nearest neighbors to consider in localerr.
#' @param eigStep The number of eigenvectors to step by during the evaluation
#' in localerr.
#' Note that this must be less than the number of samples in localerr. Default = 10.
#' @param colIdInd The ID of the column that has the analyte ID's for the
#' independent variable. If blank, then the existing ID's are used.
#' @param colIdOut The ID of the column that has the analyte ID's for the
#' outcome variable. If blank, then the existing ID's are used.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @export
DoTestSetupAndPrediction <- function(inputDataTrain, inputDataTest, model, stype,
                                     outcomeType = 1,
                                     independentVarType = 2,
                                     metaFeatureList = c("pdf","interactionpval", "interactioncoef", "analytecoef", "localerr"),
                                     k = 2, eigStep = 1,
                                     colIdInd = "databaseId",
                                     colIdOut = "databaseId",
                                     covar = c(),
                                     useCutoff = FALSE){
  # Calculate prediction cutoffs.
  predictionCutoffs <- CalculatePredictionCutoffs(modelResults = modelResults)
  
  # Run pairwise prediction for training data.
  pred.cv <- MultiOmicsGraphPrediction::RunPairwisePrediction(inputResults = model@model.input@model.properties,
                                                              inputData = inputDataTrain,
                                                              stype = stype,
                                                              covar = covar,
                                                              independentVarType = independentVarType,
                                                              outcomeType = outcomeType)
  pred.cv.test <- MultiOmicsGraphPrediction::RunPairwisePrediction(inputResults = model@model.input@model.properties,
                                                                   inputData = inputDataTest,
                                                                   stype = stype,
                                                                   covar = covar,
                                                                   independentVarType = independentVarType,
                                                                   outcomeType = outcomeType)
  
  # Get metafeatures for testing.
  metafeaturesTest <- MultiOmicsGraphPrediction::GetTestMetaFeatures(predictionsTrain = pred.cv,
                                                                     predictionsTest = pred.cv.test,
                                                                     stype = stype,
                                                                     inputDataTrain = inputDataTrain,
                                                                     inputDataTest = inputDataTest,
                                                                     metaFeatureList = metaFeatureList,
                                                                     k = k, eigStep = eigStep,
                                                                     modelStats = model@model.input@model.properties,
                                                                     colIdInd = colIdInd,
                                                                     colIdOut = colIdOut)
  
  # Predict testing values.
  predictions <- Predict(pairs = model@pairs,
                         inputData = data$testing, 
                         metafeatures = metafeaturesTest, 
                         model = model,
                                useCutoff = useCutoff,
                                minCutoff = predictionCutoffs$min,
                                maxCutoff = predictionCutoffs$max)
  return(predictions)
}

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
#' @param pairs A list of pairs to include in the composite model.
#' @param inputData The input testing data, which should be of an 
#' IntLimData class type.
#' @param metafeatures The metafeatures calculated using GetTestMetaFeatures()
#' @param model The learned model, an object of the modelResults class.
#' @param minCutoff Mininum cutoff for the prediction.
#' @param maxCutoff Maximum cutoff for the prediction.
#' @param useCutoff Whether or not to use the cutoff for prediction. Default is FALSE.
#' @param useActivation Whether or not to use the activation function if the phenotype
#' to predict is a factor. Default is TRUE. We set to FALSE during training.
#' @export
Predict <- function(pairs, inputData, metafeatures, model, minCutoff, maxCutoff, useCutoff = FALSE,
                    useActivation = TRUE){
  
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
  if(model@model.input@outcome.type == "categorical" && useActivation == TRUE){
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