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
    quantiles <- quantile(pred, c(lowerPercentileLimit, upperPercentileLimit))
    which_quantile <- intersect(which(pred > quantiles[1]),
                                which(pred < quantiles[2]))
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
ComputePDFImportance <- function(predictions, outliers){
  
  # Compute the importance for each sample.
  importanceAll <- lapply(1:nrow(predictions), function(samp){
    # Prediction
    pred <- predictions[samp,]
    
    # Remove outliers.
    non_outliers <- which(outliers[[samp]] == 0)
    
    # Compute importance based on density.
    d <- density(pred[non_outliers])
    importance[non_outliers] <- unlist(lapply(pred[non_outliers], function(pred){
      return(d$y[which.min(abs(d$x - pred))])
    }))
    
    # Scale importance.
    importance <- importance / sum(importance)
    
    # Return.
    cat(".")
    return(importance)
  })
  return(importanceAll)
}

# Function for obtaining prediction names.
GetPredictorRaMPNames <- function(preds, analyteMappingsType1){
  # Get prediction names.
  geneNames <- unlist(lapply(colnames(preds), function(pair){
    return(strsplit(pair, "__")[[1]][1])
  }))
  metabNames <- unlist(lapply(colnames(preds), function(pair){
    return(strsplit(pair, "__")[[1]][2])
  }))
  
  # Map to source ID.
  metabId <- analyteMappingsType1[metabNames, "databaseId"]
  geneId <- paste0("gene_symbol:", geneNames)
  
  # Map to RaMP ID.
  geneRampId <- unlist(lapply(geneId, function(id){
    return(source$RaMPID[which(source$SourceID == id)][1])
  }))
  metabRampId <- unlist(lapply(metabId, function(id){
    return(source$RaMPID[which(source$SourceID == id)][1])
  }))
  
  # Return data frame.
  return(data.frame(gene = geneRampId, metab = metabRampId))
}

#' Compute the importance using pathway membership. Predictor pairs that share
#' at least one pathway have an importance of 1, and all others have an 
#' importance of 0. Outliers are also set to 0.
#' @param predictions Prediction data frame, where rows are samples, and
#' columns are predictors.
#' @param analytehaspathway A mapping from analytes to pathways. This
#' mapping has two columns: RaMPID and PathID
#' @param analyteNamesOut A data frame mapping the analyte identifier for the
#' outcome analyte from the IntLIM model to an identifier that can be mapped
#' to RaMP.
#' @param analyteNamesInd A data frame mapping the analyte identifier for the
#' independent variable analyte from the IntLIM model to an identifier that can be mapped
#' to RaMP.
#' @param source A data frame that maps analyte identifiers to RaMP internal
#' identifiers.
#' @param outliers A vector of outliers for each sample. Generated using 
#' FindDistributionalOutliers().
#' @return A list of vectors (one for each sample) with an importance metric.
ComputePathwayImportance <- function(predictions, analytehaspathway, analyteNamesOut,
                                     analyteNamesInd, source,
                                     outliers){
  
  # Obtain names of predictors.
  predNames <- GetPredictorRaMPNames(predictions, analyteNamesOut,
                                     analyteNamesInd)
  
  # Measure importance for each predictor.
  importance <- unlist(lapply(1:nrow(predNames), function(pair_i){
    
    # Get gene and metabolite names.
    gene <- predNames[pair_i,"gene"]
    metab <- predNames[pair_i, "metab"]
    
    # Check whether any RaMP ID is NA.
    retval <- 0
    if(!is.na(gene) && !is.na(metab)){
      
      # Determine whether shared pathway occurs.
      genePathways <- analytehaspathway$PathID[which(analytehaspathway$RaMPID == gene)]
      metabPathways <- analytehaspathway$PathID[which(analytehaspathway$RaMPID == metab)]
      intersect <- intersect(genePathways, metabPathways)
      if(length(intersect) > 0){
        retval <- 1
      }
    }
    return(retval)
  }))
  
  # Compute the importance for each sample.
  importanceAll <- lapply(1:nrow(predictions), function(samp){
    
    # Remove zeros.
    importance[which(outliers[[samp]] == 0)] <- 0
    cat(".")
    return(importance)
  })
  return(importanceAll)
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
  if(eigStepBy > dim(type1Similarity)){
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
        d <- as.dist(dist)
        dend <- hclust(d)
        coph_i <- stats::cophenetic(dend)
        coph_cor_i <- cor(coph_i, d)
        
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
#' @param eigStep The number of eigenvectors to step by during the evaluation.
#' Note that this must be less than the number of samples. Default = 10.
#' @param alphaMin The lowest value of alpha to investigate. Default = 0.
#' @param alphaMax The highest value of alpha to investigate. Default = 1.
#' @param alphaStep The value of alpha to step by during the evalutation.
#' Default = 0.1.
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
ComputeLocalPredImportance <- function(predictions, true, outliers, k,
                                       inputData, eigStep = 10, alphaMin = 0,
                                       alphaMax = 1, alphaStep = 0.1){
  
  # Find the optimal projection for best cluster separability.
  gene_sim <- ComputeCosineSimilarity(t(inputData@analyteType1))
  metab_sim <- ComputeCosineSimilarity(t(inputData@analyteType2))
  opt <- FindOptimalSubspaceClustering(type1Similarity = gene_sim, 
                                       type2Similarity = metab_sim,
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
    trueVals <- matrix(rep(true[knn], length(pred_local)),
                       nrow = length(knn),
                       ncol = ncol(pred_local))
    
    err <- abs(trueVals - pred_local)
    medErr <- unlist(lapply(1:ncol(err), function(c){
      return(median(err[,c]))
    }))
    importance <- 1 / medErr
    
    # Remove outliers.
    outliers <- which(outliers[[samp]] == 1)
    importance[outliers] <- 0

    # Scale importance.
    importance <- importance / sum(importance)
    
    # Return.
    cat(".")
    return(importance)
  })
  return(importanceAll)
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
ComputeGlobalPredImportance <- function(predictions, true, outliers){
  
  # Compute the importance for each sample.
  importanceAll <- lapply(1:nrow(predictions), function(samp){
    # Prediction
    pred <- predictions[-samp,]
    
    # Error
    trueVals <- matrix(rep(true, length(pred)),
                       nrow = length(true),
                       ncol = ncol(pred))
    
    err <- abs(trueVals - pred)
    medErr <- unlist(lapply(1:ncol(err), function(c){
      return(median(err[,c]))
    }))
    importance <- 1 / medErr
    
    # Remove outliers.
    outliers <- which(outliers[[samp]] == 1)
    importance[outliers] <- 0
    
    # Scale importance.
    importance <- importance / sum(importance)
    
    # Return.
    cat(".")
    return(importance)
  })
  return(importanceAll)
}