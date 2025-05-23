% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/consensus_clustering.R
\name{consensus_clustering}
\alias{consensus_clustering}
\title{Perform Consensus Clustering Using K-means}
\usage{
consensus_clustering(data, k, n_runs = 10)
}
\arguments{
\item{data}{A numeric matrix where rows are observations (e.g., cells) and columns are features (e.g., genes).}

\item{k}{The number of clusters to form in k-means.}

\item{n_runs}{The number of repeated k-means runs to build the consensus matrix. Default is 10.}
}
\value{
A list containing:
\describe{
  \item{consensus_matrix}{A matrix indicating the proportion of times each pair of observations was clustered together.}
}
}
\description{
This function performs consensus clustering by repeatedly applying k-means clustering
and aggregating the clustering results to produce a consensus matrix.
}
\examples{
set.seed(1)
data <- matrix(rnorm(100*10), nrow = 100, ncol = 10)
result <- consensus_clustering(data, k = 3, n_runs = 5)
str(result$consensus_matrix)
}
