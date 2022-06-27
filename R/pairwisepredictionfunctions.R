#' Given each significant pairwise model and the input data, predict the phenotype
#' for each sample. Recall that IntLIM models take the following form, where a_i
#' and a_j are a pair of analytes.
#' a_i ~ beta0 + beta1(a_j) + beta2(phenotype) + beta3(a_j:phenotype) + beta4...n(covariates)
#' Therefore, to predict phenotype given the betas learned by IntLIM, we use the
#' following model:
#' p ~ (a_i - (beta0 + beta1(a_j) + beta4...n(covariates)) / (beta2 + beta3(a_j))
#' @param inputResults The data frame of filtered results from IntLIM
#'  model and processing results (output of ProcessResults()). All results must
#'  include learned covariate weights (i.e. must be run with save.covar.pvals = TRUE)
#' @param inputData An IntLimData object that includes the input.
#' @param stype The phenotype (outcome) to predict. This can be either a categorical
#' or numeric outcome.
#' @param covar The clinical covariates to include in the model. These should be the same
#' covariates that were included when running the IntLIM linear models.
#' @param independentVarType The independent variable type (1 or 2)
#' @param outcomeType The outcome type (1 or 2)
#' @return A data frame of predictions, where each column is a sample and each
#' row is a predictor.
#' @export
RunPairwisePrediction <- function(inputResults, inputData, stype="", covar=c(),
                                  independentVarType = 2, outcomeType = 1){
  
  # Extract data needed for model.
  covariates <- inputData@sampleMetaData
  independent.vars <- data.frame()
  dependent.vars <- data.frame()
  independent.vars <- as.data.frame(inputData@analyteType1)
  if(independentVarType == 2){
    independent.vars <- as.data.frame(inputData@analyteType2)
  }
  dependent.vars <- as.data.frame(inputData@analyteType1)
  if(outcomeType == 2){
    dependent.vars <- as.data.frame(inputData@analyteType2)
  }
  
  # Extract coefficients.
  which_start <- which(colnames(inputResults) == "rsquared")[1]
  if(is.na(which_start)){
    which_start <- which(colnames(inputResults) == "FDRadjPval")[1]
  }
  which_start <- which_start + 1
  coefficients <- inputResults[,c(1:2)]
  if(which_start <= (ncol(inputResults))){
    coefficients <- inputResults[,c(1:2, which_start:(ncol(inputResults)))]
  }
  coefficients$Analyte1 <- as.character(coefficients$Analyte1)
  coefficients$Analyte2 <- as.character(coefficients$Analyte2)
  which_interact <- which(grepl(":", colnames(coefficients)) == TRUE)
  
  # Construct matrix of coefficients.
  intercept <- matrix(coefficients[,"(Intercept)"], 
                      nrow=length(coefficients[,"(Intercept)"]), 
                      ncol=ncol(independent.vars))
  ind_var <- matrix(coefficients$a, 
                    nrow=length(coefficients$a), 
                    ncol=ncol(independent.vars))
  interact_var <- matrix(coefficients[,which_interact], 
                         nrow=length(coefficients[,which_interact]), 
                         ncol=ncol(independent.vars))
  phen_var <- matrix(coefficients$type, nrow=length(coefficients$type), 
                     ncol=ncol(independent.vars))
  
  # Compute the numerator of the prediction for each subject 
  # (sans covariate terms)
  ind_term <- independent.vars[coefficients$Analyte1,] * ind_var
  dep_term <- as.matrix(dependent.vars[coefficients$Analyte2,])
  pred_phenotype <- dep_term - (intercept + ind_term)
  rownames(dep_term) <- rownames(inputResults)
  rownames(intercept) <- rownames(inputResults)
  rownames(ind_term) <- rownames(inputResults)
  
  # If there are covariates, include the covariate terms in the prediction
  # by subtracting them.
  if(length(covar) > 0) {
    all_cov_terms<- lapply(covar, function(cov_name){
      this_covariate_term <- matrix()
      
      # We expect that the terms are one-hot-encoded now, so it should be numeric.
      this_coefficient_mat <- matrix(coefficients[,cov_name], 
                                     nrow=length(coefficients[,cov_name]), 
                                     ncol=dim(covariates)[1])
      this_covariate_mat <- t(matrix(covariates[,cov_name], 
                                     ncol=length(coefficients[,cov_name]), 
                                     nrow=dim(covariates)[1]))
      this_covariate_term <- this_coefficient_mat * this_covariate_mat
      return(this_covariate_term)
    })
    
    # Include the covariates in the numerator prediction.
    final_covariate_val <- Reduce('+', all_cov_terms)
    pred_phenotype <- pred_phenotype - final_covariate_val
  }
  
  # Calculate the denominator and divide.
  div_term <- independent.vars[coefficients$Analyte1,] * interact_var
  div_term <- div_term + phen_var
  pred_phenotype <- pred_phenotype / div_term 
  rownames(div_term) <- rownames(inputResults)
  
  # For discrete phenotypes only, round the value.
  if(!is.numeric(inputData@sampleMetaData[,stype])){
    pred_phenotype[which(pred_phenotype >= 1)]
    pred_phenotype[which(pred_phenotype >= 0)]
    pred_phenotype <- round(pred_phenotype,digits=0)
  }
  pred_phenotype = as.data.frame(pred_phenotype)
  
  # Add analyte information for each prediction.
  rownames(pred_phenotype) <- paste(coefficients$Analyte1, coefficients$Analyte2,
                                    sep = "__")
  colnames(pred_phenotype) <- rownames(inputData@sampleMetaData)
  return(t(pred_phenotype))
}

#' Given a graph and a phenotype prediction for each significant pair, generate a 
#' new graph with the phenotype predictions included as the edge weights.
#' @param predictions A matrix of predictions. Each signficant pair of analytes 
#' results in a prediction for each subject.
#' @param coRegulationGraph An igraph object. This graph is the co-regulation graph
#' generated using IntLIM analysis of analyte pairs.
#' @return A list of igraph objects, one per sample, where each node is an analyte
#' and each edge is a prediction. The weight of the edge is the predicted value
#' using the pair of analytes represented by the edge. 
#' @export
ProjectPredictionsOntoGraph <- function(predictions, coRegulationGraph){
  
  # Convert graph into data frame.
  edges <- igraph::as_data_frame(coRegulationGraph, what = "edges")
  
  # Define continuous color based on prediction.
  pal <- grDevices::colorRampPalette(c("limegreen", "purple"))

  # Modify n copies of graph, where n is the number of subjects.
  new_graphs <- lapply(1:dim(predictions)[1], function(i){

    # Modify graph weight.
    subject_graph <- data.frame(edges)
    subject_graph$weight <- predictions[i,]
    
    # Modify color.
    bin_count <- 100
    preds <- predictions[i,][which(!is.na(predictions[i,]))]
    intervals <- seq(range(preds)[1], range(preds)[2], 
                     by = (range(preds)[2] - range(preds)[1])
                     / (bin_count - 1))
    subject_color_scale <- rep(NA, length(predictions[i,]))
    subject_color_scale[which(!is.na(predictions[i,]))] <- 
      findInterval(preds, intervals)
    subject_graph$color <- pal(bin_count + 1)[subject_color_scale]

    # Modify node properties.
    node_df <- igraph::as_data_frame(coRegulationGraph, what = "vertices")
    
    # Create graph.
    final_graph = igraph::graph_from_data_frame(subject_graph, vertices = node_df)
    return(final_graph)
  })
  names(new_graphs)<-rownames(predictions)
  return(new_graphs)
}

#' Perform one-hot encoding of the sample metadata. This is useful when categorical
#' covariates are used. This one-hot encoding will be used later in predictions,
#' and will simplify the prediction process because the names of the result columns
#' and the covariates will match up.
#' @param inputData An object with the following fields:
#' @param covar The clinical covariates to include in the model. These should be the same
#' covariates that were included when running the IntLIM linear models.
#' @return A list with two elements: (1) the new set of covariates, and (2)
#' the new sample meta data, as encoded.
#' @export
OneHotEncoding <- function(inputData, covar=c()){
  # Extract the covariate values from the input data.
  covariateVals <- inputData@sampleMetaData
  
  # If the covariates are factors, then one-hot encode them. Else, add the same
  # values back into the list.
  newCovars <- c()
  for(c in covar){
    if(is.numeric(covariateVals[,c])){
      newCovars <- c(newCovars, c)
    }
    # Find all factor covariates derived from the original covariate.
    for(fac in unique(covariateVals[,c])){
      # Append new one-hot-encoded factor covariates.
      covariateVals[,paste0(c, fac)] <- 0
      covariateVals[which(covariateVals[,c] == fac),paste0(c, fac)] <- 1
      
      # Append to the list of new covariates.
      newCovars <- c(newCovars, paste0(c, fac))
    }
  }
  # Return results.
  retVal <- list(covar = newCovars, sampleMetaData = covariateVals)
  return(retVal)
}