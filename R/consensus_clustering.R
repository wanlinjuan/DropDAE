#' Perform Consensus Clustering Using K-means
#'
#' This function performs consensus clustering by repeatedly applying k-means clustering
#' and aggregating the clustering results to produce a consensus matrix.
#'
#' @param data A numeric matrix where rows are observations (e.g., cells) and columns are features (e.g., genes).
#' @param k The number of clusters to form in k-means.
#' @param n_runs The number of repeated k-means runs to build the consensus matrix. Default is 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{consensus_matrix}{A matrix indicating the proportion of times each pair of observations was clustered together.}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(100*10), nrow = 100, ncol = 10)
#' result <- consensus_clustering(data, k = 3, n_runs = 5)
#' str(result$consensus_matrix)

consensus_clustering <- function(data, k, n_runs = 10) {
  ncells <- nrow(data)
  consensus_matrix <- matrix(0, nrow = ncells, ncol = ncells)
  ARI_i <- rep(NA, n_runs)

  for (i in 1:n_runs) {
    # Perform k-means clustering
    kmeans_result <- kmeans(data, centers = k, nstart = 1)
    cluster_assignments <- kmeans_result$cluster

    # Update the consensus matrix
    for (c in 1:k) {
      indices <- which(cluster_assignments == c)
      consensus_matrix[indices, indices] <- consensus_matrix[indices, indices] + 1
    }
  }

  # Normalize the consensus matrix
  consensus_matrix <- consensus_matrix / n_runs
  # mean_ARI <- mean(ARI_i)
  return(list(consensus_matrix=consensus_matrix))
}




