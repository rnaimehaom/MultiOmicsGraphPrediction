% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositemodelfunctions.R
\name{ComputeSignificance}
\alias{ComputeSignificance}
\title{Compute the significance value for a given prediction. You may use information
gain, odds ratio, or t-statistic.}
\usage{
ComputeSignificance(pred, trueVal, pruningMethod = "error.t.test")
}
\arguments{
\item{trueVal}{The true values (predictions or outcomes) of the input data.}

\item{pruningMethod}{The method to use for pruning. Right now, only "error.t.test" is valid.}

\item{pairs}{A list of pairs to include in the composite model.}

\item{includeVarianceTest}{Scale the t-score by the f-statistic (the ratio of variances).
Only applicable when the pruning method is error.t.test. Default is FALSE.}
}
\description{
Compute the significance value for a given prediction. You may use information
gain, odds ratio, or t-statistic.
}
