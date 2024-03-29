% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositemodelfunctions.R
\name{PrunePredictors}
\alias{PrunePredictors}
\title{Given multiple composite predictors, prune the predictors that are not needed.}
\usage{
PrunePredictors(
  compositeSubgraphs,
  previousModels,
  modelResults,
  verbose = FALSE,
  makePlots = FALSE,
  pruningMethod = "error.t.test",
  tolerance = 1e-05,
  modelRetention = "stringent",
  minCutoff,
  maxCutoff,
  useCutoff = FALSE,
  weights,
  individualPerformance,
  weightCutoff
)
}
\arguments{
\item{compositeSubgraphs}{A list of pairs to include in the composite model.}

\item{previousModels}{A list of the previous models that were consolidated.}

\item{modelResults}{A ModelResults object.}

\item{verbose}{Whether or not to print out each step.}

\item{makePlots}{Whether or not to plot the pruned model at each step.}

\item{pruningMethod}{The method to use for pruning. Right now, only "error.t.test" is valid.}

\item{tolerance}{Tolerance factor when computing equality of two numeric values.}

\item{modelRetention}{Strategy for model retention. "stringent" (the default)
retains only models that improve the prediction score. "lenient" also retains models that
neither improve nor reduce the prediction score.}

\item{minCutoff}{Mininum cutoff for the prediction.}

\item{maxCutoff}{Maximum cutoff for the prediction.}

\item{useCutoff}{Whether or not to use the cutoff for prediction. Default is FALSE.}

\item{weights}{The weights for each predictor, calculated using ComputeMetaFeatureWeights()}

\item{individualPerformance}{The score (using the pruning method) for each individual component of the model.}

\item{weightCutoff}{Only consider pairs with weight above this cutoff.}

\item{previousModelsTarget}{A list of the targets in the previous models that were consolidated.}

\item{previousModelsSource}{A list of the sources in the previous models that were consolidated.}
}
\value{
A list of sets, where each set is a neighborhood of nodes.
}
\description{
Given multiple composite predictors, prune the predictors that are not needed.
}
