% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coregulationgraphfunctions.R
\name{BuildGraphDataFrame}
\alias{BuildGraphDataFrame}
\title{Construct a data frame that includes all information needed to build an igraph
object. This includes the names of the two co-regulated analytes, whether the
analytes have a positive or negative co-regulation, and the shape of the analytes.}
\usage{
BuildGraphDataFrame(inputResults)
}
\arguments{
\item{inputResults}{A data frame containing the analyte names for each model
and the coefficients, p-values, and R^2 values.}
}
\description{
Construct a data frame that includes all information needed to build an igraph
object. This includes the names of the two co-regulated analytes, whether the
analytes have a positive or negative co-regulation, and the shape of the analytes.
}
