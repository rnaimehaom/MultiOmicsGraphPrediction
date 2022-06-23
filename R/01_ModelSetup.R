#' This function performs the following tasks:
#' 1. Run IntLIM to obtain all pairwise models.
#' 2. Filter the IntLIM results to obtain a subset of models.
#' 3. Predict the phenotype of all training data given each pairwise model.
#' 4. Build the co-regulation graph using the set of pairwise models.
#' 5. Initialize model with parameters.
#' @param inputData An object with the following fields:
#' @param stype The phenotype (outcome) to predict. This can be either a categorical
#' or numeric outcome.
#' @param covar The clinical covariates to include in the model. These should be the same
#' covariates that were included when running the IntLIM linear models.
#' @param independentVarType The independent variable type (1 or 2)
#' @param outcomeType The outcome type (1 or 2)
#' @param continuous Whether or not the outcome is continuous. Default is TRUE.
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
#' @param edgeTypeList List containing one or more of the following to include
#' in the line graph:
#' - "shared.outcome.analyte"
#' - "shared.independent.analyte"
#' - "analyte.chain"
#' @param learningRate Learning rate to use during training. Default is 0.2
#' @export
DoModelSetup <- function(inputData, stype,
                         outcomeType = 1,
                                independentVarType = 2,
                                continuous = TRUE,
                                pvalcutoff = 1, rsquaredCutoff = 0, interactionCoeffPercentile = 0,
                                metaFeatureList = c("pdf","interactionpval", "interactioncoef", "analytecoef", "localerr"),
                                k = 2, eigStep = 1,
                                colIdInd = "databaseId",
                                colIdOut = "databaseId",
                         edgeTypeList = c("shared.outcome.analyte", "shared.independent.analyte"),
                                learningRate = 0.2,
                         covar = c()){
  
  # Run IntLIM.
  myres <- IntLIM::RunIntLim(inputData = inputData, 
                                stype=stype,
                                save.covar.pvals = TRUE, 
                                outcome = outcomeType, 
                                independent.var.type = independentVarType,
                                covar = covar,
                                continuous = continuous)
  
  # Process the results. If user has specified the R^2 value, then filter using this
  # pre-specified value. Otherwise, infer the R^2 value as the value of maximum density.
  myres.filt <- NULL
  if(rsquaredCutoff > 0){
    myres.filt <- IntLIM::ProcessResults(inputResults = myres, inputData = inputData, 
                                         pvalcutoff = pvalcutoff, rsquaredCutoff = rsquaredCutoff, 
                                         interactionCoeffPercentile = interactionCoeffPercentile)
  }else{
    myres.all <- IntLIM::ProcessResults(inputResults = myres, inputData = inputData, 
                                         pvalcutoff = 1, rsquaredCutoff = 0, 
                                         interactionCoeffPercentile = 0)
    myres.r2.density <- density(myres.all$rsquared)
    peak.r2 <- myres.r2.density$x[which.max(myres.r2.density$y)]
    if(peak.r2 > 1){peak.r2 <- 1}
    if(peak.r2 < 0){peak.r2 <- 0}
    myres.filt <- IntLIM::ProcessResults(inputResults = myres, inputData = inputData, 
                                         pvalcutoff = pvalcutoff, rsquaredCutoff = peak.r2, 
                                         interactionCoeffPercentile = interactionCoeffPercentile)
  }
  
  # Perform one-hot-encoding on the covariates.
  encoding <- MultiOmicsGraphPrediction::OneHotEncoding(covar = covar, inputData = inputData)
  covar = encoding$covar
  inputData@sampleMetaData <- encoding$sampleMetaData

  # Run pairwise prediction.
  pred <- MultiOmicsGraphPrediction::RunPairwisePrediction(inputResults = myres.filt, 
                                                              inputData = inputData,
                                                              stype = stype,
                                                              independentVarType = independentVarType,
                                                              outcomeType = outcomeType)
  coreg <- MultiOmicsGraphPrediction::BuildCoRegulationGraph(myres.filt)
  projectedGraph <- MultiOmicsGraphPrediction::ProjectPredictionsOntoGraph(coRegulationGraph = coreg,
                                                                              predictions = pred)
  
  # Compute meta features.
  metafeatures <- MultiOmicsGraphPrediction::GetMetaFeatures(predictions = pred,
                                                                  stype = stype,
                                                                  inputData = inputData,
                                                                  metaFeatureList = metaFeatureList,
                                                                  k = k, eigStep = eigStep,
                                                                  modelStats = myres.filt,
                                                                  colIdInd = colIdInd,
                                                                  colIdOut = colIdOut)
  
  
  # Initialize model.
  stypeClass <- "categorical"
  if(continuous == TRUE){
    stypeClass <- "numeric"
  }
  modelInput <- MultiOmicsGraphPrediction::FormatInput(predictionGraphs = projectedGraph, 
                                                          coregulationGraph = coreg,
                                                          inputData = inputData, 
                                                          stype.class = stypeClass,
                                                          stype = stype,
                                                          edgeTypeList = edgeTypeList,
                                                          metaFeatures = metafeatures,
                                                          modelProperties = myres.filt,
                                                          outcome = outcomeType,
                                                          independent.var.type = independentVarType)
  modelResults <- MultiOmicsGraphPrediction::InitializeGraphLearningModel(modelInput = modelInput, learningRate = learningRate)
  
  # Return initialized mode,
  return(modelResults)
}