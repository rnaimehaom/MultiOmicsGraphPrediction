% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plottingfunctions.R
\name{SaveGraphPredictionPlots}
\alias{SaveGraphPredictionPlots}
\title{For each sample, save a plot of phenotype predictions from each edge in the
IntLIM graph. The edge weight is the phenotype prediction
using the pair of analytes connected to the edge. Weights correspond to color
in the graph. A color bar is shown for reference, and the true phenotype is
listed at the top of the plot.}
\usage{
SaveGraphPredictionPlots(
  graphWithPredictions,
  inputData,
  stype,
  dirName,
  continuous = TRUE
)
}
\arguments{
\item{graphWithPredictions}{An igraph object. This graph is the co-regulation graph
generated using IntLIM analysis of analyte pairs. Weights correspond to phenotype
predictions.}

\item{inputData}{Named list (output of FilterData()) with gene expression,
metabolite abundances, and associated meta-data}

\item{stype}{The phenotype of interest. This should correspond to a column in the
input data.}

\item{dirName}{The name of the directory where the output images will be saved.}

\item{continuous}{A boolean indicating whether the phenotype is continuous (TRUE)
or discrete (FALSE).}
}
\description{
For each sample, save a plot of phenotype predictions from each edge in the
IntLIM graph. The edge weight is the phenotype prediction
using the pair of analytes connected to the edge. Weights correspond to color
in the graph. A color bar is shown for reference, and the true phenotype is
listed at the top of the plot.
}
