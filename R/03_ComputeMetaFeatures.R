#' Finds and returns the distributional outliers for each sample.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param lowerPercentileLimit Sample-specific percentile below which outliers
#' will be removed.
#' @param upperPercentileLimit Sample-specific percentile above which outliers
#' will be removed.
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
GetAllImportanceMetrics <- function(predictions, inputData,
                                       lowerPercentileLimit = 0.1,
                                       upperPercentileLimit = 0.9,
                                       metricList = c("pdf", "localerr", "globalerr",
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
  # Find all outliers.
  outliers <- FindDistributionalOutliers(predictions = predictions, 
                                         lowerPercentileLimit = lowerPercentileLimit,
                                         upperPercentileLimit = upperPercentileLimit)
  
  # Compute metrics.
  metrics <- list()
  if("pdf" %in% metricList){
    print("Computing PDF importance")
    metrics$pdf <- ComputePDFImportance(predictions = predictions, 
                                        outliers = outliers)
  }
  if("localerr" %in% metricList){
    print("Computing Local Error importance")
    trueVals <- inputData@sampleMetaData[,stype]
    names(trueVals) <- rownames(inputData@sampleMetaData)
    metrics$localerr <- ComputeLocalErrorImportance(predictions = predictions, 
                                                    true = trueVals,
                                                    k = k,
                                                    inputData = inputData, 
                                                    eigStep = eigStep, 
                                                    alphaMin = alphaMin,
                                                    alphaMax = alphaMax, 
                                                    alphaStep = alphaStep)
  }
  if("globalerr" %in% metricList){  
    print("Computing Global Error importance")
    metrics$globalerr <- ComputeGlobalErrorImportance(predictions = predictions, 
                                                      true = inputData@sampleMetaData[,stype])
  }
  if("pathway" %in% metricList){
    print("Computing Pathway importance")
    metrics$pathway <- ComputePathwayImportance(predictions = predictions, 
                                                inputData = inputData,
                                                colIdInd = colIdInd,
                                                colIdOut = colIdOut)
  }
  if("reaction" %in% metricList){
    print("Computing Reaction importance")
    metrics$reaction <- ComputeReactionImportance(predictions = predictions, 
                                                  inputData = inputData,
                                                  colIdInd = colIdInd,
                                                  colIdOut = colIdOut)
  }
  if("interactionpval" %in% metricList){
    print("Computing Interaction Term p-Value importance")
    metrics$interactionpval <- ComputePvalImportance(predictions = predictions, 
                                                     modelStats = modelStats)
  }
  if("interactioncoef" %in% metricList){
    print("Computing Interaction Coefficient importance")
    metrics$interactioncoef <- ComputeInteractionCoefImportance(predictions = predictions, 
                                                                modelStats = modelStats)
  }
  if("analytecoef" %in% metricList){
    print("Computing Analyte Coefficient importance")
    metrics$analytecoef <- ComputeAnalyteCoefImportance(predictions = predictions, 
                                                        modelStats = modelStats)
  }
  
  return(metrics)
}

#' Finds and returns the distributional outliers for each sample.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param lowerPercentileLimit Sample-specific percentile below which outliers
#' will be removed.
#' @param upperPercentileLimit Sample-specific percentile above which outliers
#' will be removed.
#' @return A list of vectors (one for each sample) where the outliers are masked
#' with a 1. Non-outliers are labeled as 0.
FindDistributionalOutliers <- function(predictions,
                                       lowerPercentileLimit = 0.1, 
                                       upperPercentileLimit = 0.9){
  
  # Compute the mask for each sample.
  outlier_mask <- lapply(1:nrow(predictions), function(samp){
    # Prediction
    pred <- predictions[samp,]
    
    # Find quantiles.
    quantiles <- stats::quantile(pred, c(lowerPercentileLimit, upperPercentileLimit))
    which_quantile <- union(which(pred < quantiles[1]),
                                which(pred > quantiles[2]))
    mask <- rep(0, length(pred))
    mask[which_quantile] <- 1
    return(mask)
  })
  return(outlier_mask)
}

#' Computes the importance as the density in the probability density function,
#' normalized over all densities. Outliers are set to 0.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param outliers A vector of outliers for each sample. Generated using 
#' FindDistributionalOutliers().
#' @return A list of vectors (one for each sample) with an importance metric.
#' @export
ComputePDFImportance <- function(predictions, outliers){
  
  # Compute the importance for each sample.
  importanceAll <- lapply(1:nrow(predictions), function(samp){
    # Prediction
    pred <- predictions[samp,]
    
    # Remove outliers.
    non_outliers <- which(outliers[[samp]] == 0)

    # Compute importance based on density.
    importance <- rep(0, length(outliers[[samp]]))
    d <- stats::density(pred[non_outliers])
    importance[non_outliers] <- unlist(lapply(pred[non_outliers], function(pred){
      return(d$y[which.min(abs(d$x - pred))])
    }))
    
    # Scale importance to be between 0 and 1.
    importance <- importance / max(importance)
    
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
#' importance of 0. Outliers are also set to 0.
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
#' importance of 0. Outliers are also set to 0.
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

#' Computes the importance as the median absolute error for local predictors
#' (i.e. the predictions for k nearest neighbors of each sample).
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param true Named vector of the true outcome values.
#' @param outliers A vector of outliers for each sample. Generated using 
#' FindDistributionalOutliers().
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
ComputeLocalErrorImportance <- function(predictions, true, outliers, k,
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
    importance <- 1 / medErr

    # Scale importance to be between 0 and 1
    importance <- importance / max(importance)
    
    # Return.
    cat(".")
    return(importance)
  })
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}

#' Computes the importance as the median absolute error for all predictors
#' on all samples except the sample in question.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param true Named vector of the true outcome values.
#' @param outliers A vector of outliers for each sample. Generated using 
#' FindDistributionalOutliers().
#' @return A list of vectors (one for each sample) with an importance metric.
#' @export
ComputeGlobalErrorImportance <- function(predictions, true, outliers){
  
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
    importance <- 1 / medErr
    
    # Scale importance.
    importance <- importance / sum(importance)
    
    # Return.
    cat(".")
    return(importance)
  })
  importanceFinal <- t(do.call(cbind, importanceAll))
  rownames(importanceFinal) <- rownames(predictions)
  colnames(importanceFinal) <- colnames(predictions)
  return(importanceFinal)
}