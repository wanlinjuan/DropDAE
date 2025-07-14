#' Prepare Autoencoder Inputs by Simulating Dropout and Normalization
#'
#' This function prepares training and testing datasets for a denoising autoencoder by simulating dropout
#' and applying normalization methods (log normalization or SCTransform). It creates corrupted and uncorrupted
#' versions of both training and test data.
#'
#' @param training_data A gene-by-cell count matrix used for training.
#' @param test_data A gene-by-cell count matrix used for testing.
#' @param dropout_mid Midpoint parameter for the dropout function. Default is 0.5.
#' @param dropout_shape Shape parameter for the dropout function. Default is -1.
#' @param normalization Normalization method to apply. Either \code{"log"} for log normalization or \code{"sct"} for SCTransform. Default is \code{"log"}.
#' @param seed Random seed for reproducibility. Default is 1234.
#'
#' @return A list containing four matrices:
#' \describe{
#'   \item{data_train_corrupt}{Corrupted (dropout-added) training data.}
#'   \item{data_test_corrupt}{Corrupted (dropout-added) testing data.}
#'   \item{data_train_uncorrupt}{Original uncorrupted training data.}
#'   \item{data_test_uncorrupt}{Original uncorrupted testing data.}
#' }
#' @export
#'
#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData
#' @importFrom splatter newSplatParams setParams
#' @importFrom scuttle logNormCounts
#' @importFrom Seurat CreateSeuratObject SCTransform NormalizeData FindVariableFeatures ScaleData GetAssayData
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("splatter", quietly = TRUE)) {
#'   set.seed(123)
#'   library(Seurat)
#'   library(splatter)
#'   params <- newSplatParams(batchCells = 200, nGenes = 100)
#'   sim <- splatSimulate(params, verbose = FALSE)
#'   data <- as.matrix(counts(sim))
#'   result <- prepare_autoencoder_inputs(training_data = data, test_data = data, normalization = "log")
#'   str(result)
#' }


prepare_autoencoder_inputs <- function(training_data,       # count
                                       test_data,           # count
                                       dropout_mid = 0.5,
                                       dropout_shape = -1,
                                       normalization = "sct",
                                       seed = 1234) {

  set.seed(seed)

  # Ensure matrix format
  training_data <- as.matrix(training_data)
  test_data <- as.matrix(test_data)

  # Create SingleCellExperiment objects
  sce_train <- SingleCellExperiment(assays = list(counts = training_data))
  sce_test <- SingleCellExperiment(assays = list(counts = test_data))

  sce_train = logNormCounts(sce_train)
  sce_train@assays@data@listData[["TrueCounts"]] = sce_train@assays@data@listData[["counts"]]
  sce_train@assays@data@listData$CellMeans = sce_train@assays@data@listData[["counts"]]
  colData(sce_train)$Cell = colnames(sce_train)
  rowData(sce_train)$Gene = rownames(sce_train)

  sce_test = logNormCounts(sce_test)
  sce_test@assays@data@listData[["TrueCounts"]] = sce_test@assays@data@listData[["counts"]]
  sce_test@assays@data@listData$CellMeans = sce_test@assays@data@listData[["counts"]]
  colData(sce_test)$Cell = colnames(sce_test)
  rowData(sce_test)$Gene = rownames(sce_test)

  # Simulate dropout on training
  params_train <- newSplatParams(seed = seed)
  params_train <- setParams(params_train, nGenes = nrow(training_data),
                            update = list(dropout.type = "experiment",
                                          dropout.mid = dropout_mid,
                                          dropout.shape = dropout_shape,
                                          seed = seed))
  sce_train_corrupt <- splatter:::splatSimDropout(sce_train, params_train)

  # Simulate dropout on test
  # params_test <- newSplatParams(seed = seed + 1)
  # params_test <- setParams(params_test, nGenes = nrow(test_data),
  #                          update = list(dropout.type = "experiment",
  #                                        dropout.mid = dropout_mid,
  #                                        dropout.shape = dropout_shape,
  #                                        seed = seed + 1))
  # sim_test <- splatter:::splatSimDropout(sce_test, params_test)
  sce_test_corrupt <- splatter:::splatSimDropout(sce_test, params_train)


  # normalization
  seurat_train_corrupt = as.Seurat(sce_train_corrupt)
  seurat_test_corrupt = as.Seurat(sce_test_corrupt)
  seurat_train_uncorrupt = CreateSeuratObject(training_data,  assay="originalexp")
  seurat_test_uncorrupt = CreateSeuratObject(test_data, assay="originalexp")

  if(normalization == "log"){
    seurat_train_corrupt = NormalizeData(seurat_train_corrupt)
    seurat_train_corrupt = FindVariableFeatures(seurat_train_corrupt)
    seurat_train_corrupt = ScaleData(seurat_train_corrupt)

    seurat_test_corrupt = NormalizeData(seurat_test_corrupt)
    seurat_test_corrupt = FindVariableFeatures(seurat_test_corrupt)
    seurat_test_corrupt = ScaleData(seurat_test_corrupt)

    seurat_train_uncorrupt = NormalizeData(seurat_train_uncorrupt)
    seurat_train_uncorrupt = FindVariableFeatures(seurat_train_uncorrupt)
    seurat_train_uncorrupt = ScaleData(seurat_train_uncorrupt)

    seurat_test_uncorrupt = NormalizeData(seurat_test_uncorrupt)
    seurat_test_uncorrupt = FindVariableFeatures(seurat_test_uncorrupt)
    seurat_test_uncorrupt = ScaleData(seurat_test_uncorrupt)

  } else if(normalization == "sct"){
    seurat_train_corrupt = SCTransform(seurat_train_corrupt,assay="originalexp")

    seurat_test_corrupt = SCTransform(seurat_test_corrupt,assay="originalexp")

    seurat_train_uncorrupt = SCTransform(seurat_train_uncorrupt,assay="originalexp")

    seurat_test_uncorrupt = SCTransform(seurat_test_uncorrupt,assay="originalexp")
  }




  # Extract corrupted matrices
  data_train_corrupt <- GetAssayData(seurat_train_corrupt, slot = "scale.data")   # uncorrupted training
  data_test_corrupt <- GetAssayData(seurat_test_corrupt, slot = "scale.data")   # corrupted training
  data_train_uncorrupt <- GetAssayData(seurat_train_uncorrupt, slot = "scale.data")    # uncorrupted testing
  data_test_uncorrupt <- GetAssayData(seurat_test_uncorrupt, slot = "scale.data")    # corrupted testing

  return(list(data_train_corrupt = t(data_train_corrupt),
              data_test_corrupt = t(data_test_corrupt),
              data_train_uncorrupt = t(data_train_uncorrupt),
              data_test_uncorrupt = t(data_test_uncorrupt)))
}
