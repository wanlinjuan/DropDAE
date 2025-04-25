#' Generate Triplet Indices for Contrastive Learning
#'
#' This function generates triplet indices for contrastive loss training,
#' based on a consensus matrix of cell-cell similarities.
#'
#' @param consensus_matrix A numeric matrix representing cell-to-cell similarities (values between 0 and 1).
#' @param batch_indices A vector of integers indicating the indices of the current batch.
#' @param pos_threshold A threshold above which cells are considered positive pairs. Default is 0.8.
#' @param neg_threshold A threshold below which cells are considered negative pairs. Default is 0.2.
#'
#' @return A list containing:
#' \describe{
#'   \item{pos_indices}{Vector of positive indices corresponding to each anchor cell.}
#'   \item{neg_indices}{Vector of negative indices corresponding to each anchor cell.}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' consensus <- matrix(runif(100), nrow = 10)
#' batch <- 1:10
#' triplets <- generate_triplet_indices(consensus, batch)
#' str(triplets)

generate_triplet_indices <- function(consensus_matrix, batch_indices, pos_threshold = 0.8, neg_threshold = 0.2) {
  pos_indices <- sapply(batch_indices, function(i) {
    high_similarity <- which(consensus_matrix[i, ] > pos_threshold)
    valid_indices <- intersect(high_similarity, batch_indices)
    if (length(valid_indices) > 1) sample(valid_indices[valid_indices != i], 1) else i
  })

  neg_indices <- sapply(batch_indices, function(i) {
    low_similarity <- which(consensus_matrix[i, ] < neg_threshold)
    valid_indices <- intersect(low_similarity, batch_indices)
    if (length(valid_indices) > 0) sample(valid_indices, 1) else i
  })

  return(list(pos_indices = pos_indices, neg_indices = neg_indices))
}
