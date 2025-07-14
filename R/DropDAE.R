#' DropDAE: Denoising Autoencoder with Triplet Loss for scRNA-seq Data
#'
#' This function trains a denoising autoencoder using corrupted scRNA-seq data,
#' with optional triplet loss based on consensus clustering to improve clustering and imputation.
#'
#' @param training_data Matrix of raw training counts.
#' @param test_data Matrix of raw testing counts.
#' @param normalization Normalization method: "sct" or "log". Default is "sct".
#' @param num_epochs Number of training epochs. Default is 300.
#' @param batch_size Mini-batch size for training. Default is 1000.
#' @param lr Initial learning rate. Default is 0.01.
#' @param hidden_size Vector specifying hidden layer sizes. Default is c(64, 32, 64).
#' @param lambda Weighting for the triplet loss component. Default is 10.
#' @param k Number of clusters used for consensus clustering. Default is 2.
#' @param patience Number of epochs with no improvement before early stopping. Default is 5.
#' @param margin Margin for triplet loss. Default is 10.
#' @param num_initializations Number of random initializations to try. Default is 1.
#' @param apply_early_stopping Whether to apply early stopping. Default is TRUE.
#' @param lr_scheduler Whether to apply learning rate decay. Default is TRUE.
#' @param select_genes Optional vector of genes to subset before training.
#' @param activation Activation function: "tanh" or "sigmoid" or "relu". Default is "tanh".
#' @param weight_decay Whether to apply L2 regularization. Default is TRUE.
#' @param cc_runs Number of consensus clustering runs. Default is 50.
#' @param BN Whether to apply batch normalization. Default is TRUE.
#' @param gradient_clip Whether to apply gradient clipping. Default is TRUE.
#' @param clip_value Maximum allowed gradient norm if clipping is used. Default is 5.0.
#'
#' @return A list containing:
#' \describe{
#'   \item{best_x_hat}{Matrix of denoised and imputed test data.}
#'   \item{loss_data_all_df}{Data frame recording loss curves across training epochs.}
#' }
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes ggtitle xlab ylab theme
#' @importFrom stats kmeans
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scuttle logNormCounts
#' @importFrom splatter newSplatParams setParams
#' @importFrom torch torch_float32 nn_linear nn_batch_norm1d nn_identity nn_relu torch_tensor
#' @importFrom torch nn_sigmoid nn_tanh nn_mse_loss torch_norm torch_relu
#' @importFrom torch tensor_dataset nn_utils_clip_grad_norm_ nn_sequential optim_adam dataloader
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom Seurat ScaleData SCTransform GetAssayData
#'
#' @examples
#' \dontrun{
#' # Example using simulated scRNA-seq data
#' }

DropDAE = function(training_data,       # count
                   test_data,           # count
                   normalization="sct",
                   num_epochs=300,
                   batch_size=1000,
                   lr=0.01,
                   hidden_size=c(64,32,64),
                   lambda=10,
                   k=2,
                   patience=5,
                   margin=10,
                   num_initializations=1,
                   # init_weights = "xavier",           # default or xavier
                   apply_early_stopping=TRUE,
                   lr_scheduler=TRUE,                   # autotune learning rate
                   select_genes=NULL,
                   activation = "relu",                 # sigmoid, tanh
                   weight_decay = TRUE,                 # L2 regularization
                   cc_runs = 50,                        # number of clusterings running for consensus clustering
                   BN = TRUE,                          # batch normalization
                   gradient_clip = TRUE,    # Whether to apply gradient clipping
                   clip_value = 5.0,         # Maximum allowed norm for gradients
                   seed=seed
) {

  set.seed(seed)

  best_loss <- Inf  # Track the best loss
  best_x_hat <- NULL  # Track the best reconstructed output
  loss_data_all <- list()  # To store loss curves for each initialization

  input_all = prepare_autoencoder_inputs(training_data,       # count
                                         test_data,
                                         normalization=normalization)
  tr1 = input_all$data_train_uncorrupt
  tr2 = input_all$data_train_corrupt
  te1 = input_all$data_test_uncorrupt
  te2 = input_all$data_test_corrupt

  common_i1 = intersect(colnames(tr1),colnames(tr2))
  common_i2 = intersect(colnames(te1),colnames(te2))
  common_i = intersect(common_i1,common_i2)

  tr1 = tr1_full = tr1[,common_i]
  tr2 = tr2_full = tr2[,common_i]
  te1 = te1_full = te1[,common_i]
  te2 = te2_full = te2[,common_i]

  # Subset data for selected genes, if provided
  if (!is.null(select_genes)) {
    cat("Subsetting data for selected genes...\n")

    # Assuming column names of tr1/tr2/te1/te2 are the gene names
    gene_indices <- which(common_i %in% select_genes)

    if (length(gene_indices) == 0) {
      stop("None of the specified genes are found in the data.")
    }

    # Subset the data for the selected genes
    tr1 <- tr1[, gene_indices]
    tr2 <- tr2[, gene_indices]
    te1 <- te1[, gene_indices]
    te2 <- te2[, gene_indices]

    # Update the number of genes (features)
    p <- length(gene_indices)
    cat(sprintf("Number of selected genes: %d\n", p))
  } else {
    cat("Processing all genes...\n")
    p <- dim(tr1)[2]  # Number of features (genes)
  }

  # Convert input data to torch tensors
  tr1 <- torch_tensor(as.matrix(tr1), dtype = torch_float32())
  tr2 <- torch_tensor(as.matrix(tr2), dtype = torch_float32())
  te1 <- torch_tensor(as.matrix(te1), dtype = torch_float32())
  te2 <- torch_tensor(as.matrix(te2), dtype = torch_float32())

  # Get the dimensions of the data
  ncells <- dim(tr1)[1]
  p <- dim(tr1)[2]
  cell_indices = torch_tensor(1:ncells)

  # Scheduler parameters
  initial_lr <- lr
  step_size <- 50  # Decay the learning rate every 50 epochs
  gamma <- 0.1     # Multiply learning rate by 0.1 at each step

  # if(init_weights == "xavier"){
  #   initialize_weights <- function(m) {
  #     if (inherits(m, "nn_linear")) {
  #       nn_init_xavier_uniform_(m$weight)
  #       if (!is.null(m$bias)) {
  #         nn_init_constant_(m$bias, 0)
  #       }
  #     }
  #   }
  # }
  #

  for (init in 1:num_initializations) {
    cat(sprintf("Initialization %d/%d\n", init, num_initializations))

    # Reset best loss and early stopping variables for this initialization
    best_loss <- Inf
    epochs_no_improve <- 0
    losses <- rep(NA, num_epochs)  # Store combined loss for each epoch
    mse_losses <- rep(NA, num_epochs)  # Store MSE loss for each epoch
    triplet_losses <- rep(NA, num_epochs)  # Store triplet loss for each epoch

    # Define the autoencoder with custom initialization
    encoder <- nn_sequential(
      nn_linear(p, hidden_size[1]),
      if (BN) nn_batch_norm1d(hidden_size[1]) else nn_identity(),
      if (activation == "relu") nn_relu() else if (activation == "sigmoid") nn_sigmoid() else nn_tanh(),
      nn_linear(hidden_size[1], hidden_size[2]),
      if (BN) nn_batch_norm1d(hidden_size[2]) else nn_identity(),
      if (activation == "relu") nn_relu() else if (activation == "sigmoid") nn_sigmoid() else nn_tanh()
    )

    decoder <- nn_sequential(
      nn_linear(hidden_size[2], hidden_size[1]),
      if (BN) nn_batch_norm1d(hidden_size[1]) else nn_identity(),
      if (activation == "relu") nn_relu() else if (activation == "sigmoid") nn_sigmoid() else nn_tanh(),
      nn_linear(hidden_size[1], p)
    )


    mse_loss <- nn_mse_loss(reduction = "mean")
    triplet_loss <- function(anchors, positives, negatives, margin) {
      pos_dist <- torch_norm(anchors - positives, p=2)
      neg_dist <- torch_norm(anchors - negatives, p=2)
      loss <- torch_relu(pos_dist - neg_dist + margin)
      return(loss)
    }

    autoencoder <- nn_sequential(encoder, decoder)
    # autoencoder$apply(initialize_weights)
    if(weight_decay){
      optimizer <- optim_adam(autoencoder$parameters, lr = lr, weight_decay = 0.0001)
    }else{
      optimizer <- optim_adam(autoencoder$parameters, lr = lr)
    }


    # Create dataloader
    dataset <- tensor_dataset(tr2, tr1, torch_tensor(1:ncells))
    dataloader <- dataloader(dataset, batch_size = batch_size, shuffle = TRUE)


    # Training loop
    for (epoch in 1:num_epochs) {

      # Adjust learning rate at the start of each step
      if (lr_scheduler && epoch %% step_size == 0) {
        new_lr <- initial_lr * gamma^(epoch %/% step_size)
        optimizer$param_groups[[1]]$lr <- new_lr
        initial_lr <- new_lr  # Update initial_lr
        cat(sprintf("Epoch %d: Adjusted learning rate to %.6f\n", epoch, new_lr))
      }


      batch_losses <- c()
      batch_mse_losses <- c()
      batch_triplet_losses <- c()
      autoencoder$train()

      if (lambda > 0) {
        # Perform clustering using updated bottleneck representations
        bottleneck_encodings <- encoder(tr2)$detach()
        consensus_run <- consensus_clustering(as.matrix(bottleneck_encodings), k = k, n_runs = cc_runs)
        consensus_matrix <- consensus_run$consensus_matrix
      }

      # Iterate over batches
      dataloader_iterator <- dataloader$.iter()
      for (i in 1:(length(dataloader) * batch_size / batch_size)) {
        # Retrieve batch
        batch <- tryCatch({
          dataloader_iterator$.next()
        }, error = function(e) {
          dataloader_iterator <- dataloader$.iter()
          dataloader_iterator$.next()
        })

        input <- batch[[1]]  # Corrupted data
        target <- batch[[2]]  # Uncorrupted data
        batch_indices <- as.array(batch[[3]])  # Batch cell indices

        optimizer$zero_grad()

        # Forward pass
        output <- autoencoder$forward(input)
        mse_loss_val <- mse_loss(output, target)
        batch_mse_losses <- c(batch_mse_losses, mse_loss_val$item())

        # Skip triplet loss calculation if lambda = 0
        if (lambda > 0) {

          # if(cc_runs>1){
          # Generate triplet indices using the consensus matrix
          triplet_indices <- generate_triplet_indices(consensus_matrix, batch_indices, pos_threshold = 0.8, neg_threshold = 0.2)
          pos_indices <- triplet_indices$pos_indices
          neg_indices <- triplet_indices$neg_indices

          # Compute triplet loss
          anchors <- bottleneck_encodings[batch_indices, ]
          positives <- bottleneck_encodings[pos_indices, ]
          negatives <- bottleneck_encodings[neg_indices, ]
          triplet_loss_val <- triplet_loss(anchors, positives, negatives, margin)
          batch_triplet_losses <- c(batch_triplet_losses, triplet_loss_val$item())
        } else {
          triplet_loss_val <- torch_tensor(0)  # No triplet loss
        }

        # Combined loss
        total_loss <- mse_loss_val + lambda * triplet_loss_val
        total_loss$backward()
        if (gradient_clip) {
          nn_utils_clip_grad_norm_(autoencoder$parameters, max_norm = clip_value)
        }
        optimizer$step()

        batch_losses <- c(batch_losses, total_loss$item())
      }

      # Track epoch loss
      epoch_loss <- mean(batch_losses)
      epoch_mse_loss <- mean(batch_mse_losses)
      epoch_triplet_loss <- mean(batch_triplet_losses)

      losses[epoch] <- epoch_loss
      mse_losses[epoch] <- epoch_mse_loss
      if (lambda > 0) triplet_losses[epoch] <- epoch_triplet_loss

      # Evaluate on the test set
      autoencoder$eval()
      x_hat <- autoencoder$forward(te2)
      x_hat <- as.data.frame(as.matrix(x_hat))

      # Check for the best loss for this initialization
      if (epoch_loss < best_loss) {
        best_loss <- epoch_loss
        best_x_hat <- x_hat
        epochs_no_improve <- 0
      } else {
        epochs_no_improve <- epochs_no_improve + 1
      }

      # Apply early stopping if enabled
      if (apply_early_stopping && epochs_no_improve >= patience) {
        cat(sprintf("Early stopping at epoch %d for initialization %d\n", epoch, init))
        break
      }
      cat(sprintf("Epoch %d/%d, Loss: %.6f, MSE Loss: %.6f, Triplet Loss: %.6f\n",
                  epoch, num_epochs, epoch_loss, epoch_mse_loss, epoch_triplet_loss))
    }

    # Store loss data for this initialization
    loss_data_all[[init]] <- data.frame(
      Epoch = 1:num_epochs,
      Loss = losses,
      MSE_Loss = mse_losses,
      Triplet_Loss = triplet_losses,
      Initialization = paste("Init", init)
    )
  }


  # Filter the data for early stopping
  loss_data_all_df <- do.call(rbind, loss_data_all) %>%
    dplyr::filter(!is.na(Loss)) %>%
    dplyr::group_by(Initialization) %>%
    dplyr::mutate(MaxEpoch = max(Epoch)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Epoch <= MaxEpoch)

  # Plot filtered loss curves
  loss_plot = ggplot(loss_data_all_df, aes(x = Epoch, y = Loss, color = Initialization)) +
    geom_line(na.rm = TRUE) +
    theme_minimal() +
    ggtitle("Filtered Loss Curves for Multiple Initializations") +
    xlab("Epoch") +
    ylab("Loss") +
    theme(legend.position = "right")
  print(loss_plot)


  # Plot filtered mse_loss curves
  mse_loss_plot = ggplot(loss_data_all_df, aes(x = Epoch, y = MSE_Loss, color = Initialization)) +
    geom_line(na.rm = TRUE) +
    theme_minimal() +
    ggtitle("Filtered mse_loss Curves for Multiple Initializations") +
    xlab("Epoch") +
    ylab("mse_loss") +
    theme(legend.position = "right")
  print(mse_loss_plot)


  if (lambda > 0) {
    # Plot filtered triplet_loss curves
    triplet_loss_plot = ggplot(loss_data_all_df, aes(x = Epoch, y = Triplet_Loss, color = Initialization)) +
      geom_line(na.rm = TRUE) +
      theme_minimal() +
      ggtitle("Filtered triplet_loss Curves for Multiple Initializations") +
      xlab("Epoch") +
      ylab("triplet_loss") +
      theme(legend.position = "right")
    print(triplet_loss_plot)
  }

  cat(sprintf("Best initialization achieved a loss of %.6f\n", best_loss))


  if (!is.null(select_genes)) {
    cat("Adjusting best_x_hat back to full input size...\n")
    full_x_hat <- as.matrix(te1_full)
    best_x_hat <- as.matrix(best_x_hat)
    full_x_hat[, gene_indices] <- best_x_hat
    best_x_hat <- full_x_hat
  }

  return(list(best_x_hat = best_x_hat, loss_data_all_df = loss_data_all_df))
}
