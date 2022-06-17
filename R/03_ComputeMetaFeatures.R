#' Computes the metafeatures for each sample and model.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param metricList A list of the valid metrics to include. Valid metrics are
#' "pdf", "localerr", "globalerr", and "pathway".
#' @param k The number of nearest neighbors to consider in localerr.
#' @param inputData The input data read in using the function IntLIM function
#' ReadData.
#' @param eigStep The number of eigenvectors to step by during the evaluation
#' in localerr.
#' Note that this must be less than the number of samples in localerr. Default = 10.
#' @param alphaMin The lowest value of alpha to investigate in localerr. Default = 0.
#' @param alphaMax The highest value of alpha to investigate in localerr. Default = 1.
#' @param alphaStep The value of alpha to step by during the evaluation in localerr.
#' Default = 0.1.
#' @param analyteNamesOut A data frame mapping the analyte identifier for the
#' outcome analyte from the IntLIM model to an identifier that can be mapped
#' to RaMP. Used in pathway importance.
#' @param analyteNamesInd A data frame mapping the analyte identifier for the
#' independent variable analyte from the IntLIM model to an identifier that can be mapped
#' to RaMP. Used in pathway importance.
#' @param modelStats A data frame that includes the interaction p-values and
#' interaction coefficients for each pair (such as the one output by IntLIM's
#' ProcessResults function)
#' @param stype Phenotype or outcome to use in models.
#' @param colIdInd The ID of the column that has the analyte ID's for the
#' independent variable. If blank, then the existing ID's are used.
#' @param colIdOut The ID of the column that has the analyte ID's for the
#' outcome variable. If blank, then the existing ID's are used.
#' @return A list of data frames (one for each sample) with predictor importance
#' measured according to the listed criteria (one column per metric, one row
#' per predictor).
#' @export
GetMetaFeatures <- function(predictions, inputData,
                                       lowerPercentileLimit = 0.1,
                                       upperPercentileLimit = 0.9,
                                       metaFeatureList = c("pdf", "localerr", "globalerr",
                                                      "pathway", "reaction",
                                                      "interactionpval", "interactioncoef",
                                                      "analytecoef"),
                                       modelStats = "",
                                       k = k,
                                       eigStep = 10, 
                                       alphaMin = 0,
                                       alphaMax = 1, 
                                       alphaStep = 0.1,
                                       analyteNamesOut = "",
                                       analyteNamesInd = "",
                                    stype = "",
                                    colIdInd = "",
                                    colIdOut = ""){
  
  # Compute metrics.
  metaFeatures <- list()
  if("pdf" %in% metaFeatureList){
    print("Computing PDF importance")
    metaFeatures$pdf <- ComputePDFImportance(predictions = predictions)
  }
  if("localerr" %in% metaFeatureList){
    print("Computing Local Error importance")
    trueVals <- inputData@sampleMetaData[,stype]
    names(trueVals) <- rownames(inputData@sampleMetaData)
    metaFeatures$localerr <- ComputeLocalErrorImportance(predictions = predictions, 
                                                    true = trueVals,
                                                    k = k,
                                                    inputData = inputData, 
                                                    eigStep = eigStep, 
                                                    alphaMin = alphaMin,
                                                    alphaMax = alphaMax, 
                                                    alphaStep = alphaStep)
  }
  if("globalerr" %in% metaFeatureList){  
    print("Computing Global Error importance")
    metaFeatures$globalerr <- ComputeGlobalErrorImportance(predictions = predictions, 
                                                      true = inputData@sampleMetaData[,stype])
  }
  if("pathway" %in% metaFeatureList){
    print("Computing Pathway importance")
    metaFeatures$pathway <- ComputePathwayImportance(predictions = predictions, 
                                                inputData = inputData,
                                                colIdInd = colIdInd,
                                                colIdOut = colIdOut)
  }
  if("reaction" %in% metaFeatureList){
    print("Computing Reaction importance")
    metaFeatures$reaction <- ComputeReactionImportance(predictions = predictions, 
                                                  inputData = inputData,
                                                  colIdInd = colIdInd,
                                                  colIdOut = colIdOut)
  }
  if("interactionpval" %in% metaFeatureList){
    print("Computing Interaction Term p-Value importance")
    metaFeatures$interactionpval <- ComputePvalImportance(predictions = predictions, 
                                                     modelStats = modelStats)
  }
  if("interactioncoef" %in% metaFeatureList){
    print("Computing Interaction Coefficient importance")
    metaFeatures$interactioncoef <- ComputeInteractionCoefImportance(predictions = predictions, 
                                                                modelStats = modelStats)
  }
  if("analytecoef" %in% metaFeatureList){
    print("Computing Analyte Coefficient importance")
    metaFeatures$analytecoef <- ComputeAnalyteCoefImportance(predictions = predictions, 
                                                        modelStats = modelStats)
  }
  
  return(metaFeatures)
}

#' Computes the metafeatures for each sample and model in the testing data.
#' @param predictionsTrain Prediction data frame for training data, where rows are samples, and
#' columns are predictors.
#' @param predictionsTest Prediction data frame for testing data, where rows are samples, and
#' columns are predictors.
#' @param metricList A list of the valid metrics to include. Valid metrics are
#' "pdf", "localerr", "globalerr", and "pathway".
#' @param k The number of nearest neighbors to consider in localerr.
#' @param inputDataTrain The training data (in IntLimData format).
#' @param inputDataTest The testing data (in IntLimData format).
#' @param eigStep The number of eigenvectors to step by during the evaluation
#' in localerr.
#' Note that this must be less than the number of samples in localerr. Default = 10.
#' @param alphaMin The lowest value of alpha to investigate in localerr. Default = 0.
#' @param alphaMax The highest value of alpha to investigate in localerr. Default = 1.
#' @param alphaStep The value of alpha to step by during the evaluation in localerr.
#' Default = 0.1.
#' @param analyteNamesOut A data frame mapping the analyte identifier for the
#' outcome analyte from the IntLIM model to an identifier that can be mapped
#' to RaMP. Used in pathway importance.
#' @param analyteNamesInd A data frame mapping the analyte identifier for the
#' independent variable analyte from the IntLIM model to an identifier that can be mapped
#' to RaMP. Used in pathway importance.
#' @param modelStats A data frame that includes the interaction p-values and
#' interaction coefficients for each pair (such as the one output by IntLIM's
#' ProcessResults function)
#' @param stype Phenotype or outcome to use in models.
#' @param colIdInd The ID of the column that has the analyte ID's for the
#' independent variable. If blank, then the existing ID's are used.
#' @param colIdOut The ID of the column that has the analyte ID's for the
#' outcome variable. If blank, then the existing ID's are used.
#' @return A list of data frames (one for each sample) with predictor importance
#' measured according to the listed criteria (one column per metric, one row
#' per predictor).
#' @export
GetTestMetaFeatures <- function(predictionsTrain, predictionsTest, inputDataTrain, inputDataTest,
                                    lowerPercentileLimit = 0.1,
                                    upperPercentileLimit = 0.9,
                                metaFeatureList = c("pdf", "localerr", "globalerr",
                                                   "pathway", "reaction",
                                                   "interactionpval", "interactioncoef",
                                                   "analytecoef"),
                                    modelStats = "",
                                    k = k,
                                    eigStep = 10, 
                                    alphaMin = 0,
                                    alphaMax = 1, 
                                    alphaStep = 0.1,
                                    analyteNamesOut = "",
                                    analyteNamesInd = "",
                                    stype = "",
                                    colIdInd = "",
                                    colIdOut = ""){
  
  # Compute metrics.
  metaFeatures <- list()
  if("pdf" %in% metaFeatureList){
    print("Computing PDF importance")
    metaFeatures$pdf <- ComputePDFImportance(predictions = predictionsTest)
  }
  if("localerr" %in% metaFeatureList){
    print("Computing Local Error importance")
    trueVals <- inputData@sampleMetaData[,stype]
    names(trueVals) <- rownames(inputData@sampleMetaData)
    metaFeatures$localerr <- ComputeLocalErrorImportanceTest(predictions = predictionsTrain, 
                                                    true = trueVals,
                                                    k = k,
                                                    inputDataTrain = inputDataTrain, 
                                                    inputDataTest = inputDataTest,
                                                    eigStep = eigStep, 
                                                    alphaMin = alphaMin,
                                                    alphaMax = alphaMax, 
                                                    alphaStep = alphaStep)
  }
  if("globalerr" %in% metaFeatureList){  
    print("Computing Global Error importance")
    metaFeatures$globalerr <- ComputeGlobalErrorImportance(predictions = predictionsTrain, 
                                                      true = inputData@sampleMetaData[,stype])
    # Adjust dimensionality for test data.
    metaFeatures$globalerr <- t(matrix(rep(metaFeatures$globalerr[1,], nrow(predictionsTest),
                                         ncol = nrow(predictionsTest))))
    rownames(metaFeatures$globalerr) <- rownames(predictionsTest)
    colnames(metaFeatures$globalerr) <- colnames(predictionsTest)
  }
  if("pathway" %in% metaFeatureList){
    print("Computing Pathway importance")
    metaFeatures$pathway <- ComputePathwayImportance(predictions = predictionsTest, 
                                                inputData = inputDataTest,
                                                colIdInd = colIdInd,
                                                colIdOut = colIdOut)
    # Adjust dimensionality for test data.
    metaFeatures$pathway <- t(matrix(rep(metaFeatures$pathway[1,], nrow(predictionsTest),
                                                 ncol = nrow(predictionsTest))))
    rownames(metaFeatures$pathway) <- rownames(predictionsTest)
    colnames(metaFeatures$pathway) <- colnames(predictionsTest)
  }
  if("reaction" %in% metaFeatureList){
    print("Computing Reaction importance")
    metaFeatures$reaction <- ComputeReactionImportance(predictions = predictionsTest, 
                                                  inputData = inputDataTest,
                                                  colIdInd = colIdInd,
                                                  colIdOut = colIdOut)
    # Adjust dimensionality for test data.
    metaFeatures$reaction <- t(matrix(rep(metaFeatures$reaction[1,], nrow(predictionsTest),
                                           ncol = nrow(predictionsTest))))
    rownames(metaFeatures$reaction) <- rownames(predictionsTest)
    colnames(metaFeatures$reaction) <- colnames(predictionsTest)
  }
  if("interactionpval" %in% metaFeatureList){
    print("Computing Interaction Term p-Value importance")
    metaFeatures$interactionpval <- ComputePvalImportance(predictions = predictionsTrain, 
                                                     modelStats = modelStats)
    metaFeatures$interactionpval <- t(matrix(rep(metaFeatures$interactionpval[1,], nrow(predictionsTest),
                                       ncol = nrow(predictionsTest))))
    rownames(metaFeatures$interactionpval) <- rownames(predictionsTest)
    colnames(metaFeatures$interactionpval) <- colnames(predictionsTest)
  }
  if("interactioncoef" %in% metaFeatureList){
    print("Computing Interaction Coefficient importance")
    metaFeatures$interactioncoef <- ComputeInteractionCoefImportance(predictions = predictionsTrain, 
                                                                modelStats = modelStats)
    # Adjust dimensionality for test data.
    metaFeatures$interactioncoef <- t(matrix(rep(metaFeatures$interactioncoef[1,], nrow(predictionsTest),
                                          ncol = nrow(predictionsTest))))
    rownames(metaFeatures$interactioncoef) <- rownames(predictionsTest)
    colnames(metaFeatures$interactioncoef) <- colnames(predictionsTest)
  }
  if("analytecoef" %in% metaFeatureList){
    print("Computing Analyte Coefficient importance")
    metaFeatures$analytecoef <- ComputeAnalyteCoefImportance(predictions = predictionsTrain, 
                                                        modelStats = modelStats)
    
    # Adjust dimensionality for test data.
    metaFeatures$analytecoef <- t(matrix(rep(metaFeatures$analytecoef[1,], nrow(predictionsTest),
                                                 ncol = nrow(predictionsTest))))
    rownames(metaFeatures$analytecoef) <- rownames(predictionsTest)
    colnames(metaFeatures$analytecoef) <- colnames(predictionsTest)
  }
  
  return(metaFeatures)
}

#' Computes the importance as the density in the probability density function,
#' normalized over all densities.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @return A list of vectors (one for each sample) with an importance metric.
#' @export
ComputePDFImportance <- function(predictions){
  
  # Compute the importance for each sample.
  importanceAll <- lapply(1:nrow(predictions), function(samp){
    # Prediction
    pred <- predictions[samp,]

    # Compute importance based on density.
    importance <- rep(0, length(pred))
    d <- stats::density(pred)
    importance <- unlist(lapply(pred, function(pred){
      return(d$y[which.min(abs(d$x - pred))])
    }))
    
    # Get the density percentile over all pairs.
    percentileVals <- seq(1, 100, by = 1) / 100
    quantiles <- stats::quantile(importance, percentileVals)
    percentiles <- unlist(lapply(importance, function(c){
      cvec <- rep(c, length(quantiles))
      which_match <- which.min(abs(cvec - quantiles))
      return(percentileVals[which_match])
    }))
    importance <- percentiles
    
    # Return.
    cat(".")
    return(importance)
  })
  importance <- t(do.call(cbind, importanceAll))
  rownames(importance) <- rownames(predictions)
  colnames(importance) <- colnames(predictions)
  return(importance)
}

#' Function for obtaining database identifiers for each predictor.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param inputData Input data, which may include database ID mappings for
#' the analytes.
#' @param colIdInd The ID of the column that has the analyte ID's for the
#' independent variable. If blank, then the existing ID's are used.
#' @param colIdOut The ID of the column that has the analyte ID's for the
#' outcome variable. If blank, then the existing ID's are used.
GetPredictorIDs <- function(preds, inputData, colIdInd, colIdOut){
  # Get prediction names.
  indAnalytes <- unlist(lapply(colnames(preds), function(pair){
    return(strsplit(pair, "__")[[1]][1])
  }))
  outAnalytes <- unlist(lapply(colnames(preds), function(pair){
    return(strsplit(pair, "__")[[1]][2])
  }))
  
  # Map to source ID.
  indId <- indAnalytes
  outId <- outAnalytes
  if(length(inputData@analyteType1MetaData) > 0 && colIdInd != "" && indAnalytes[[1]] %in% rownames(inputData@analyteType1MetaData)){
    indId <- inputData@analyteType1MetaData[indAnalytes, colIdInd]
  }else if(length(inputData@analyteType2MetaData) > 0 && colIdInd != "" && indAnalytes[[1]] %in% rownames(inputData@analyteType2MetaData)){
    indId <- inputData@analyteType2MetaData[indAnalytes, colIdInd]
  }
  if(length(inputData@analyteType1MetaData) > 0 && colIdOut != "" && outAnalytes[[1]] %in% rownames(inputData@analyteType1MetaData)){
    outId <- inputData@analyteType1MetaData[outAnalytes, colIdOut]
  }else if(length(inputData@analyteType2MetaData) > 0 && colIdOut != "" && outAnalytes[[1]] %in% rownames(inputData@analyteType2MetaData)){
    outId <- inputData@analyteType2MetaData[outAnalytes, colIdOut]
  }
  
  # Return data frame.
  return(data.frame(indId = indId, outId = outId))
}

#' Compute the importance using pathway membership. Predictor pairs that share
#' at least one pathway have an importance of 1, and all others have an 
#' importance of 0.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param inputData Input data, which may include database ID mappings for
#' the analytes.
#' @param colIdInd The ID of the column that has the analyte ID's for the
#' independent variable. If blank, then the existing ID's are used.
#' @param colIdOut The ID of the column that has the analyte ID's for the
#' outcome variable. If blank, then the existing ID's are used.
#' @return A list of vectors (one for each sample) with an importance metric.
ComputePathwayImportance <- function(predictions, inputData, colIdInd,
                                     colIdOut){
  
  # Obtain names of predictors.
  predNames <- GetPredictorIDs(predictions, inputData, colIdInd, colIdOut)

  # Find which pairs share pathways.
  sharesPathway <- unlist(lapply(1:nrow(predNames), function(i){
    shares <- FALSE
    pwayResult <- RaMP::getPathwayFromAnalyte(c(predNames[i,"indId"], predNames[i,"outId"]))
    pwayResultInd <- pwayResult$pathwayId[which(pwayResult$inputId == predNames[i,"indId"])]
    pwayResultOut <- pwayResult$pathwayId[which(pwayResult$inputId == predNames[i,"outId"])]
    if(length(intersect(pwayResultInd, pwayResultOut)) > 0){
      shares <- TRUE
    }
    return(shares)
  }))
  importance <- rep(0, length(sharesPathway))
  importance[which(sharesPathway == TRUE)] <- 1
  
  # Make copies for each sample.
  importanceAll <- lapply(rownames(predictions), function(samp){return(importance)})
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Compute the importance using reaction membership. Predictor pairs that share
#' at least one reaction have an importance of 1, and all others have an 
#' importance of 0.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param inputData Input data, which may include database ID mappings for
#' the analytes.
#' @param colIdInd The ID of the column that has the analyte ID's for the
#' independent variable. If blank, then the existing ID's are used.
#' @param colIdOut The ID of the column that has the analyte ID's for the
#' outcome variable. If blank, then the existing ID's are used.
#' @return A list of vectors (one for each sample) with an importance metric.
ComputeReactionImportance <- function(predictions, inputData, colIdInd,
                                      colIdOut){
  # Obtain names of predictors.
  predNames <- GetPredictorIDs(predictions, inputData, colIdInd, colIdOut)
  
  # Find which pairs share reactions.
  sharesRxn <- unlist(lapply(1:nrow(predNames), function(i){
    shares <- FALSE
    tryCatch({
      rxnResult <- rampFastCata(predNames[i,"outId"])$rxn_partner_ids
      rxnResultAll <- unlist(lapply(rxnResult, function(res){
        return(strsplit(res, "; ")[[1]])
      }))
      if(predNames[i,"indId"] %in% rxnResultAll){
        shares <- TRUE
      }
    }, error = function(cond){})
    return(shares)
  }))
  importance <- rep(0, length(sharesRxn))
  importance[which(sharesRxn == TRUE)] <- 1
  
  # Make copies for each sample.
  importanceAll <- lapply(rownames(predictions), function(samp){return(importance)})
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Compute the importance using p-value of model interaction term. Specifically,
#' the p-value percentile is calculated, and the weight is taken as 1 - the percentile,
#' so that the lowest p-values are weighted more heavily.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param modelStats A data frame that includes the interaction p-values and
#' interaction coefficients for each pair (such as the one output by IntLIM's
#' ProcessResults function)
#' @return A list of vectors (one for each sample) with an importance metric.
ComputePvalImportance <- function(predictions, modelStats){
  
  # Measure importance for each predictor.
  pvals <- unlist(lapply(colnames(predictions), function(predictor){
    
    # Get gene and metabolite names.
    srcAnalyte <- strsplit(predictor, "__")[[1]][1]
    tgtAnalyte <- strsplit(predictor, "__")[[1]][2]
    
    # Get p-value.
    return(modelStats[intersect(which(modelStats$Analyte1 == srcAnalyte), 
                                which(modelStats$Analyte2 == tgtAnalyte)),
                      "FDRadjPval"])
  }))
  
  # Get the p-value percentile over all pairs.
  percentileVals <- seq(1, 100, by = 1) / 100
  quantiles <- stats::quantile(pvals, percentileVals)
  pvalPercentiles <- unlist(lapply(pvals, function(p){
    pvec <- rep(p, length(quantiles))
    which_match <- which.min(abs(pvec - quantiles))
    return(percentileVals[which_match])
  }))
  importance <- 1 - pvalPercentiles
  
  # Make copies for each sample.
  importanceAll <- lapply(rownames(predictions), function(samp){return(importance)})
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Compute the importance using absolute value of interaction term coefficient. Specifically,
#' the coefficient percentile is calculated and taken as the weight.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param modelStats A data frame that includes the interaction p-values and
#' interaction coefficients for each pair (such as the one output by IntLIM's
#' ProcessResults function)
#' @return A list of vectors (one for each sample) with an importance metric.
ComputeInteractionCoefImportance <- function(predictions, modelStats){
  
  # Measure importance for each predictor.
  coefs <- unlist(lapply(colnames(predictions), function(predictor){
    
    # Get gene and metabolite names.
    srcAnalyte <- strsplit(predictor, "__")[[1]][1]
    tgtAnalyte <- strsplit(predictor, "__")[[1]][2]
    
    # Get coefficient.
    return(modelStats[intersect(which(modelStats$Analyte1 == srcAnalyte), 
                                which(modelStats$Analyte2 == tgtAnalyte)),
                      "interaction_coeff"])
  }))
  
  # Get the coefficient percentile over all pairs.
  percentileVals <- seq(1, 100, by = 1) / 100
  quantiles <- stats::quantile(abs(coefs), percentileVals)
  coefPercentiles <- unlist(lapply(abs(coefs), function(c){
    cvec <- rep(c, length(quantiles))
    which_match <- which.min(abs(cvec - quantiles))
    return(percentileVals[which_match])
  }))
  importance <- coefPercentiles
  
  # Make copies for each sample.
  importanceAll <- lapply(rownames(predictions), function(samp){return(importance)})
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Compute the importance using absolute value of analyte term coefficient. Specifically,
#' the coefficient percentile is calculated and taken as the weight.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param modelStats A data frame that includes the interaction p-values,
#' interaction coefficients, and other coefficients for each pair (such as the one output by IntLIM's
#' ProcessResults function)
#' @return A list of vectors (one for each sample) with an importance metric.
ComputeAnalyteCoefImportance <- function(predictions, modelStats){
  
  # Measure importance for each predictor.
  coefs <- unlist(lapply(colnames(predictions), function(predictor){
    
    # Get gene and metabolite names.
    srcAnalyte <- strsplit(predictor, "__")[[1]][1]
    tgtAnalyte <- strsplit(predictor, "__")[[1]][2]
    
    # Get coefficient.
    return(modelStats[intersect(which(modelStats$Analyte1 == srcAnalyte), 
                                which(modelStats$Analyte2 == tgtAnalyte)),
                      "type"])
  }))
  
  # Get the coefficient percentile over all pairs.
  percentileVals <- seq(1, 100, by = 1) / 100
  quantiles <- stats::quantile(abs(coefs), percentileVals)
  coefPercentiles <- unlist(lapply(abs(coefs), function(c){
    cvec <- rep(c, length(quantiles))
    which_match <- which.min(abs(cvec - quantiles))
    return(percentileVals[which_match])
  }))
  importance <- coefPercentiles
  
  # Make copies for each sample.
  importanceAll <- lapply(rownames(predictions), function(samp){return(importance)})
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Finds the optimal subspace clustering (i.e. using Grassmann Manifold Technique
#' - see PMID 30329022) given two different modalities of data (e.g. gene and metabolite).
#' It optimizes the cophenetic correlation of the hierarchicial clustering
#' of the samples using a grid search.
#' @param type1Similarity A cosine similarity matrix for the first data type, 
#' found using ComputeCosineSimilarity.
#' @param type2Similarity A cosine similarity matrix for the second data type,
#' found using ComputeCosineSimilarity.
#' @param eigStep The number of eigenvectors to step by during the evaluation.
#' Note that this must be less than the number of samples. Default = 10.
#' @param alphaMin The lowest value of alpha to investigate. Default = 0.
#' @param alphaMax The highest value of alpha to investigate. Default = 1.
#' @param alphaStep The value of alpha to step by during the evalutation.
#' Default = 0.1.
#' @return A named list including the data projected onto the merged subspace,
#' the optimal number of eigenvectors, the optimal alpha value, the clustering
#' coefficient, and the dendrogram.
FindOptimalSubspaceClustering <- function(type1Similarity, type2Similarity,
                                          eigStep = 10, alphaMin = 0,
                                          alphaMax = 1, alphaStep = 0.1){
  # Check input parameters.
  if(eigStep > nrow(type1Similarity)){
    stop("Must set eigStep to a value lower than the number of samples!")
  }
  else if(alphaMin < 0 || alphaMin > 1 || alphaMax < 0 || alphaMax > 1 ){
    stop("Both alphaMin and alphaMax must be between 0 and 1!")
  }
  else if(alphaStep > alphaMax - alphaMin){
    stop(paste("alphaStep is too big! There must be at least one step between",
               "alphaMin and alphaMax"))
  }
  else{
    # Initialize.
    L_mod <- NULL
    dendro <- NULL
    optimal_eigens <- NULL
    optimal_alpha <- NULL
    coph_cor <- NULL
    
    # Do grid search.
    for(eig in seq(from=1,to=dim(type1Similarity)[1],by=eigStep)){
      for(a in seq(from=alphaMin,to=alphaMax,by=alphaStep)){
        
        # Find clustering with the current grid values.
        L_mod_i <- CombineSubspaces(type1Similarity = type1Similarity, 
                                    type2Similarity = type2Similarity, 
                                    eigenCount = eig, alpha = a)
        
        # Compute cophenetic correlation.
        rownames(L_mod_i) <- rownames(type1Similarity)
        sim <- ComputeCosineSimilarity(L_mod_i)
        dist <- max(sim) - sim
        rownames(dist) <- rownames(sim)
        d <- stats::as.dist(dist)
        dend <- stats::hclust(d)
        coph_i <- stats::cophenetic(dend)
        coph_cor_i <- stats::cor(coph_i, d)
        
        # If cophenetic correlation is higher than previous, save.
        if(is.null(coph_cor) || coph_cor_i > coph_cor){
          coph_cor <- coph_cor_i
          L_mod <- L_mod_i
          dendro <- dend
          optimal_eigens <- eig
          optimal_alpha <- a
        }
      }
    }
    return(list("L_mod"=L_mod, "eigenvector_ct"=optimal_eigens, "alpha"=optimal_alpha,
                "cluster_coeff"=coph_cor, "dendrogram"=dendro))
  }
}

#' Combine subspaces (i.e. using Grassmann Manifold Technique
#' - see PMID 30329022) given two different modalities of data (e.g. gene and metabolite)
#' and the alpha value and the desired number of eigenvectors.
#' @param type1Similarity A cosine similarity matrix for the first data type, 
#' found using ComputeCosineSimilarity.
#' @param type2Similarity A cosine similarity matrix for the second data type,
#' found using ComputeCosineSimilarity.
#' @param eigenCount The number of eigenvectors to use.
#' @param alpha The value of alpha to use.
#' @return A named list including the data projected onto the merged subspace,
#' the optimal number of eigenvectors, the optimal alpha value, the clustering
#' coefficient, and the dendrogram.
CombineSubspaces <- function(type1Similarity, type2Similarity, eigenCount, alpha){
  
  # Make graphs.
  D_dat1 <- diag(colSums(type1Similarity))
  D_dat2 <- diag(colSums(type2Similarity))
  
  # Get normalized Laplacian for each graph.
  L_dat1 <- as.matrix(D_dat1 - type1Similarity)
  L_dat2 <- as.matrix(D_dat2 - type2Similarity)
  
  # Normalize.
  D_neg_half_dat1 <- diag(diag(1 / sqrt(abs(D_dat1))))
  D_neg_half_dat2 <- diag(diag(1 / sqrt(abs(D_dat2))))
  L_dat1 <- D_neg_half_dat1 %*% L_dat1 %*% D_neg_half_dat1
  L_dat2 <- D_neg_half_dat2 %*% L_dat2 %*% D_neg_half_dat2
  
  # Get eigenvectors for each Laplacian.
  U_dat1 <- eigen(L_dat1)$vectors[,c(1:eigenCount)]
  U_dat2 <- eigen(L_dat2)$vectors[,c(1:eigenCount)]
  
  # Combine
  L_mod <- L_dat1 + L_dat2 - alpha * ((U_dat1 %*% t(U_dat1)) + (U_dat2 %*% t(U_dat2)))
  return(L_mod)
}

#' Combine subspaces (i.e. using Grassmann Manifold Technique
#' - see PMID 30329022) given two different modalities of data (e.g. gene and metabolite)
#' and the alpha value and the desired number of eigenvectors.
#' @param type1SimilarityTrain A cosine similarity matrix for the first data type, 
#' found using ComputeCosineSimilarity on the training data.
#' @param type2SimilarityTrain A cosine similarity matrix for the second data type,
#' found using ComputeCosineSimilarity on the training data.
#' @param type1SimilarityTest A cosine similarity matrix for the first data type, 
#' found using ComputeCosineSimilarity on the testing data.
#' @param type2SimilarityTest A cosine similarity matrix for the second data type,
#' found using ComputeCosineSimilarity on the testing data.
#' @param subspaceTraining The subspace 
#' @param eigenCount The number of eigenvectors to use.
#' @param alpha The value of alpha to use.
#' @return A named list including the data projected onto the merged subspace,
#' the optimal number of eigenvectors, the optimal alpha value, the clustering
#' coefficient, and the dendrogram.
CombineSubspacesTest <- function(type1SimilarityTrain, type2SimilarityTrain, 
                                 type1SimilarityTest, type2SimilarityTest, 
                                 eigenCount, alpha){
  
  # Make graphs. We compute the diagonals on the training data only.
  D_dat1 <- diag(colSums(type1SimilarityTrain))
  D_dat2 <- diag(colSums(type2SimilarityTrain))
  L_dat1_train <- as.matrix(D_dat1 - type1SimilarityTrain)
  L_dat2_train <- as.matrix(D_dat2 - type2SimilarityTrain)
  D_neg_half_dat1 <- diag(diag(1 / sqrt(abs(D_dat1))))
  D_neg_half_dat2 <- diag(diag(1 / sqrt(abs(D_dat2))))
  L_dat1_train <- D_neg_half_dat1 %*% L_dat1_train %*% D_neg_half_dat1
  L_dat2_train <- D_neg_half_dat2 %*% L_dat2_train %*% D_neg_half_dat2
  
  # Get eigenvectors for each Laplacian on the training data.
  U_dat1_train <- eigen(L_dat1_train)$vectors[,c(1:eigenCount)]
  U_dat2_train <- eigen(L_dat2_train)$vectors[,c(1:eigenCount)]
  
  # Do not subtract from D_dat for testing data because none of the
  # data points are on the diagonal. Simply negate the values.
  L_dat1 <- as.matrix(0 - type1SimilarityTest)
  L_dat2 <- as.matrix(0 - type2SimilarityTest)
  
  # Normalize. Here, transposition is necessary to obtain the correct
  # dimensionality of the matrix.
  D_neg_half_dat1 <- diag(diag(1 / sqrt(abs(D_dat1))))
  D_neg_half_dat2 <- diag(diag(1 / sqrt(abs(D_dat2))))
  L_dat1 <- t(t(D_neg_half_dat1 %*% t(L_dat1)) %*% D_neg_half_dat1)
  L_dat2 <- t(t(D_neg_half_dat2 %*% t(L_dat2)) %*% D_neg_half_dat2)
  
  # Project testing data onto eigenvectors.
  U_dat1 <- t(L_dat1) %*% U_dat1_train
  U_dat2 <- t(L_dat2) %*% U_dat2_train
  
  # Combine
  L_mod <- L_dat1 + L_dat2 - alpha * (t(U_dat1 %*% t(U_dat1_train)) + t(U_dat2 %*% t(U_dat2_train)))
  return(L_mod)
}

#' Compute cosine similarity between samples on an adjacency matrix.
#' @param R The adjacency matrix.
#' @return A matrix of sample rows and sample columns, with the cosine
#' similarities listed as the values.
ComputeCosineSimilarity <- function(R){
  # Normalize matrix by the norm of its columns.
  euclidean_norm <- sqrt(rowSums(R^2))
  euclidean_norm_mat <- matrix(rep(euclidean_norm, ncol(R)), ncol=ncol(R))
  R_norm <- as.matrix(R / euclidean_norm_mat)
  
  # Compute matrix transpose to get cosine similarity.
  sim <- R_norm %*% t(R_norm)
  rownames(sim) <- rownames(R)
  colnames(sim) <- rownames(R)
  return(sim)
}

#' Compute cosine similarity between a testing sample and each training sample
#' in the combined subspace matrix.
#' @param opt The optimal model obtained using FindOptimalSubspaceClustering on
#' the training data.
#' @param inputDataTest The testing data.
#' @param inputDataTrain The training data.
#' @param samp The name of the testing sample against which to compute similarity.
#' @return A matrix of similarities, where the rows are the training data samples
#' and the columns are the testing data samples.
CosineSimilarityToTrainingOnSubspace <- function(opt, inputDataTest, inputDataTrain,
                                                 samp){
  
  # Compute similarity between the training samples.
  type1SimilarityTrain <- ComputeCosineSimilarity(t(inputDataTrain@analyteType1))
  type2SimilarityTrain <- ComputeCosineSimilarity(t(inputDataTrain@analyteType2))
  
  # Compute adjacency of the sample to each of the training samples using
  # the parameters of the optimal model.
  type1SimilarityTest <- CosineSimilarityToTrainingOnSingleDataType(trainingData = inputDataTrain@analyteType1,
                                                           testingData = inputDataTest@analyteType1,
                                                           samp = samp)
  type2SimilarityTest <- CosineSimilarityToTrainingOnSingleDataType(trainingData = inputDataTrain@analyteType2,
                                                         testingData = inputDataTest@analyteType2,
                                                         samp = samp)

  R <- CombineSubspacesTest(type1SimilarityTrain = type1SimilarityTrain,
                            type2SimilarityTrain = type2SimilarityTrain,
                            type1SimilarityTest = t(type1SimilarityTest), 
                        type2SimilarityTest = t(type2SimilarityTest), 
                        eigenCount = opt$eigenvector_ct, 
                        alpha = opt$alpha)
  colnames(R) <- samp
  rownames(R) <- rownames(inputDataTrain@sampleMetaData)
  
  # Compute cosine similarity on the new subspace.
  similarity <- CosineSimilarityToTrainingOnSingleDataType(trainingData = opt$L_mod,
                                                           testingData = R,
                                                           samp = samp)
  colnames(similarity) <- samp
  rownames(similarity) <- rownames(inputDataTrain@sampleMetaData)

  return(similarity)
}

#' Compute cosine similarity between a testing sample and each training sample
#' for a single data type.
#' @param opt The optimal model obtained using FindOptimalSubspaceClustering on
#' the training data.
#' @param testingData The testing data.
#' @param trainingData The training data.
#' @param samp The name of the testing sample against which to compute similarity.
#' @return A matrix of similarities, where the rows are the training data samples
#' and the columns are the testing data samples.
CosineSimilarityToTrainingOnSingleDataType <- function(trainingData, testingData,
                                                       samp){
  # Compute cosine similarity to each training sample.
  similarityList <- matrix(unlist(lapply(1:ncol(trainingData), function(c){
    dotProduct <- sum(trainingData[,c] * testingData[,samp])
    normTrain <- sqrt(sum(trainingData[,c]^2))
    normTest <- sqrt(sum(testingData[,samp]^2))
    return(dotProduct / (normTrain * normTest))
  })))

  # Return similarity as a matrix.
  return(similarityList)
}

#' Computes the importance as the median absolute error for local predictors
#' (i.e. the predictions for k nearest neighbors of each sample).
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param true Named vector of the true outcome values.
#' @param k The number of nearest neighbors to consider.
#' @param inputData The input data read in using the function IntLIM function
#' ReadData.
#' @param eigStep The number of eigenvectors to step by during the evaluation.
#' Note that this must be less than the number of samples. Default = 10.
#' @param alphaMin The lowest value of alpha to investigate. Default = 0.
#' @param alphaMax The highest value of alpha to investigate. Default = 1.
#' @param alphaStep The value of alpha to step by during the evalutation.
#' Default = 0.1.
#' @return A list of vectors (one for each sample) with an importance metric.
#' @export
ComputeLocalErrorImportance <- function(predictions, true, k,
                                       inputData, eigStep = 10, alphaMin = 0,
                                       alphaMax = 1, alphaStep = 0.1){
  
  # Find the optimal projection for best cluster separability.
  type1sim <- ComputeCosineSimilarity(t(inputData@analyteType1))
  type2sim <- ComputeCosineSimilarity(t(inputData@analyteType2))
  opt <- FindOptimalSubspaceClustering(type1Similarity = type1sim, 
                                       type2Similarity = type2sim,
                                       eigStep = eigStep, alphaMin = alphaMin,
                                       alphaMax = alphaMax, alphaStep = alphaStep)
  
  # Compute the importance for each sample.
  importanceAll <- lapply(1:nrow(predictions), function(samp){
    # Prediction
    pred <- predictions[-samp,]
    
    # Find nearest neighbors.
    similarities <- ComputeCosineSimilarity(opt$L_mod)[samp,]
    sorted_sims <- names(similarities)[order(-similarities)]
    pred_sorted_sims <- sorted_sims[which(sorted_sims %in% rownames(pred))]
    knn <- pred_sorted_sims[1:k]
    pred_local <- pred[knn,]
    
    # Error
    trueVals <- matrix(rep(true[knn], ncol(pred_local)),
                       nrow = length(knn),
                       ncol = ncol(pred_local))

    err <- abs(trueVals - pred_local)
    medErr <- unlist(lapply(1:ncol(err), function(c){
      return(stats::median(err[,c]))
    }))
    importance <- 1 - medErr

    # Get the 1 - error percentile over all pairs.
    percentileVals <- seq(1, 100, by = 1) / 100
    quantiles <- stats::quantile(importance, percentileVals)
    percentiles <- unlist(lapply(importance, function(c){
      cvec <- rep(c, length(quantiles))
      which_match <- which.min(abs(cvec - quantiles))
      return(percentileVals[which_match])
    }))
    importance <- percentiles
    
    # Return.
    cat(".")
    return(importance)
  })
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Computes the importance as the median absolute error for local predictors
#' (i.e. the predictions for k nearest neighbors of each sample).
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param true Named vector of the true outcome values.
#' @param k The number of nearest neighbors to consider.
#' @param inputDataTrain The training data (in IntLimData format).
#' @param inputDataTest The testing data (in IntLimData format).
#' @param eigStep The number of eigenvectors to step by during the evaluation.
#' Note that this must be less than the number of samples. Default = 10.
#' @param alphaMin The lowest value of alpha to investigate. Default = 0.
#' @param alphaMax The highest value of alpha to investigate. Default = 1.
#' @param alphaStep The value of alpha to step by during the evalutation.
#' Default = 0.1.
#' @return A list of vectors (one for each sample) with an importance metric.
#' @export
ComputeLocalErrorImportanceTest <- function(predictions, true, k,
                                        inputDataTrain, inputDataTest, 
                                        eigStep = 10, alphaMin = 0,
                                        alphaMax = 1, alphaStep = 0.1){
  
  # Find the optimal projection for best cluster separability on the training data.
  type1sim <- ComputeCosineSimilarity(t(inputDataTrain@analyteType1))
  type2sim <- ComputeCosineSimilarity(t(inputDataTrain@analyteType2))
  opt <- FindOptimalSubspaceClustering(type1Similarity = type1sim, 
                                       type2Similarity = type2sim,
                                       eigStep = eigStep, alphaMin = alphaMin,
                                       alphaMax = alphaMax, alphaStep = alphaStep)
  
  # Compute the importance for each sample in the testing data.
  importanceAll <- lapply(rownames(inputDataTest@sampleMetaData), function(samp){

    # Training data predictions
    pred <- predictions
    
    # Find nearest neighbors.
    similarities <- CosineSimilarityToTrainingOnSubspace(opt = opt, 
                                                         inputDataTrain = inputDataTrain, 
                                                         inputDataTest = inputDataTest, 
                                                         samp = samp)
    sorted_sims <- rownames(similarities)[order(-similarities)]
    pred_sorted_sims <- sorted_sims[which(sorted_sims %in% rownames(pred))]
    knn <- pred_sorted_sims[1:k]
    pred_local <- pred[knn,]
    
    # Error
    trueVals <- matrix(rep(true[knn], ncol(pred_local)),
                       nrow = length(knn),
                       ncol = ncol(pred_local))
    
    err <- abs(trueVals - pred_local)
    medErr <- unlist(lapply(1:ncol(err), function(c){
      return(stats::median(err[,c]))
    }))
    importance <- 1 - medErr
    
    # Get the 1 - error percentile over all pairs.
    percentileVals <- seq(1, 100, by = 1) / 100
    quantiles <- stats::quantile(importance, percentileVals)
    percentiles <- unlist(lapply(importance, function(c){
      cvec <- rep(c, length(quantiles))
      which_match <- which.min(abs(cvec - quantiles))
      return(percentileVals[which_match])
    }))
    importance <- percentiles
    
    # Return.
    cat(".")
    return(importance)
  })
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(inputDataTest@sampleMetaData)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Computes the importance as the median absolute error for all predictors
#' on all samples except the sample in question.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param true Named vector of the true outcome values.
#' @return A list of vectors (one for each sample) with an importance metric.
#' @export
ComputeGlobalErrorImportance <- function(predictions, true){
  
  # Compute the importance for each sample.
  importanceAll <- lapply(1:nrow(predictions), function(samp){
    # Prediction
    pred <- predictions[-samp,]
    
    # Error
    trueVals <- matrix(rep(true, length(pred)),
                       nrow = length(true),
                       ncol = ncol(pred))
    trueVals <- trueVals[-samp,]
    
    err <- abs(trueVals - pred)
    medErr <- unlist(lapply(1:ncol(err), function(c){
      return(stats::median(err[,c]))
    }))
    importance <- 1 - medErr
    
    # Get the 1 - error percentile over all pairs.
    percentileVals <- seq(1, 100, by = 1) / 100
    quantiles <- stats::quantile(importance, percentileVals)
    percentiles <- unlist(lapply(importance, function(c){
      cvec <- rep(c, length(quantiles))
      which_match <- which.min(abs(cvec - quantiles))
      return(percentileVals[which_match])
    }))
    importance <- percentiles
    
    # Return.
    cat(".")
    return(importance)
  })
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}