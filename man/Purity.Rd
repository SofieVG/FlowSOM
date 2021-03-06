% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_buildSOM.R
\name{Purity}
\alias{Purity}
\title{Calculate mean weighted cluster purity}
\usage{
Purity(realClusters, predictedClusters, weighted = TRUE)
}
\arguments{
\item{realClusters}{array with real cluster values}

\item{predictedClusters}{array with predicted cluster values}

\item{weighted}{logical. Should the mean be weighted
depending on the number of poins in the predicted 
clusters}
}
\value{
Mean purity score, worst score, number of clusters with score < 0.75
}
\description{
Calculate mean weighted cluster purity
}
\examples{
# Generate some random data as an example
realClusters <- sample(1:5, 100, replace = TRUE)
predictedClusters <- sample(1:6, 100, replace = TRUE)

# Calculate the FMeasure
Purity(realClusters, predictedClusters)
}
