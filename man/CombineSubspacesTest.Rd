% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metafeaturefunctions.R
\name{CombineSubspacesTest}
\alias{CombineSubspacesTest}
\title{Combine subspaces (i.e. using Grassmann Manifold Technique
- see PMID 30329022) given two different modalities of data (e.g. gene and metabolite)
and the alpha value and the desired number of eigenvectors.}
\usage{
CombineSubspacesTest(
  type1SimilarityTrain,
  type2SimilarityTrain,
  type1SimilarityTest,
  type2SimilarityTest,
  eigenCount,
  alpha
)
}
\arguments{
\item{type1SimilarityTrain}{A cosine similarity matrix for the first data type, 
found using ComputeCosineSimilarity on the training data.}

\item{type2SimilarityTrain}{A cosine similarity matrix for the second data type,
found using ComputeCosineSimilarity on the training data.}

\item{type1SimilarityTest}{A cosine similarity matrix for the first data type, 
found using ComputeCosineSimilarity on the testing data.}

\item{type2SimilarityTest}{A cosine similarity matrix for the second data type,
found using ComputeCosineSimilarity on the testing data.}

\item{eigenCount}{The number of eigenvectors to use.}

\item{alpha}{The value of alpha to use.}

\item{subspaceTraining}{The subspace}
}
\value{
A named list including the data projected onto the merged subspace,
the optimal number of eigenvectors, the optimal alpha value, the clustering
coefficient, and the dendrogram.
}
\description{
Combine subspaces (i.e. using Grassmann Manifold Technique
- see PMID 30329022) given two different modalities of data (e.g. gene and metabolite)
and the alpha value and the desired number of eigenvectors.
}
