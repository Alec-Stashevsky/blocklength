# NPPI Algorithm for Optimal Block Length in Block Bootstrap

# Function to estimate optimal block length using NPPI algorithm
nppi_optimal_block_length <- function(data, stat_function, r = 1, a = 1, initial_block_size = NULL, num_bootstrap = 1000, c_1 = 1L) {
  n <- length(data)

  # Step 1: Set initial block size if not provided
  if (is.null(initial_block_size)) {
    initial_block_size <- floor(c_1*n^(1 / (r + 4)))
  }

  # Set c_2 to 1 if r=1 and to 0.1 otherwise
  c_2 <- ifelse(r == 1, 1, 0.1)

  # Function to create overlapping blocks
  create_overlapping_blocks <- function(data, block_size) {
    N <- length(data) - block_size + 1
    blocks <- embed(data, block_size)[N:1, ]
    return(blocks)
  }

  # Step 1: Manual Moving Block Bootstrap (MBB) Resampling
  mbb_resample <- function(blocks, num_samples) {
    N <- nrow(blocks)
    resampled_blocks <- replicate(num_samples, blocks[sample(1:N, N, replace = TRUE), ], simplify = FALSE)
    return(resampled_blocks)
  }

  # Step 2: Bias Estimation
  bias_estimator <- function(data, block_size) {
    blocks_l1 <- create_overlapping_blocks(data, block_size)
    blocks_2l1 <- create_overlapping_blocks(data, 2 * block_size)

    resampled_l1 <- mbb_resample(blocks_l1, num_bootstrap)
    resampled_2l1 <- mbb_resample(blocks_2l1, num_bootstrap)

    theta_l1 <- sapply(resampled_l1, function(blocks) stat_function(as.vector(t(blocks))))
    theta_2l1 <- sapply(resampled_2l1, function(blocks) stat_function(as.vector(t(blocks))))

    # Bias = 2 * (phi_n(l1) - phi_n(2l1))
    phi_hat_n_of_l <- mean(theta_l1)
    phi_hat_n_of_2l <- mean(theta_2l1)
    bias <- 2 * (phi_hat_n_of_l - phi_hat_n_of_2l)

    return(list(bias = bias, resampled_blocks = resampled_l1, original_blocks = blocks_l1, phi_hat_n = phi_hat_n_of_l))
  }

  # Step 3: Variance Estimation using Specialized JAB
  jab_variance_estimator <- function(resampled_blocks, original_blocks, stat_function, m, phi_hat_n) {
    N <- nrow(original_blocks)
    K <- length(resampled_blocks)
    M <- N - m + 1

    # Initialize JAB point values
    jab_point_values <- numeric(M)

    for (i in 1:M) {
      blocks_to_delete <- original_blocks[i:(i + m - 1), ]
      I_star_indices <- c()

      # Identify non-overlapping bootstrap samples
      for (k in 1:K) {
        bootstrap_sample <- resampled_blocks[[k]]

        overlaps <- any(apply(blocks_to_delete, 1, function(del_block) {
          any(apply(bootstrap_sample, 1, function(sample_block) all(sample_block == del_block)))
        }))

        if (!overlaps) {
          I_star_indices <- c(I_star_indices, k)
        }
      }

      # Compute the i-th JAB point value using non-overlapping bootstrap samples
      if (length(I_star_indices) > 0) {
        non_overlapping_samples <- resampled_blocks[I_star_indices]
        jab_point_values[i] <- stat_function(as.vector(t(do.call(rbind, non_overlapping_samples))))
      } else {
        jab_point_values[i] <- NA  # Handle cases with no non-overlapping samples
      }
    }

    # Remove NA values and calculate JAB variance
    jab_point_values <- na.omit(jab_point_values)
    jab_var <- (m / (N - m)) * (1 / M) * sum((jab_point_values - phi_hat_n)^2)

    return(jab_var)
  }

  # Step 4: Compute Bias and Variance
  bias_result <- bias_estimator(data, initial_block_size)
  bias <- bias_result$bias
  resampled_blocks <- bias_result$resampled_blocks
  original_blocks <- bias_result$original_blocks
  phi_hat_n <- bias_result$phi_hat_n  # Reuse phi_hat_n from bias estimation

  m <- floor(c_2 * n^(1/3) * initial_block_size^(2/3))  # Suggested m from paper
  variance <- jab_variance_estimator(resampled_blocks, original_blocks, stat_function, m, phi_hat_n)

  # Step 5: Calculate Optimal Block Length
  C_hat_1 <- n * l^(-r) * variance
  C_hat_2 <- l * bias

  optimal_block_length <- ((2 * C_hat_2^2) / (r * C_hat_1)^(1 / (r + 2)) * n^(1 / (r + 2)))

  return(list(optimal_block_length = optimal_block_length, bias = bias, variance = variance))

}







# Example usage with a simple mean function
# stat_function <- function(x) mean(x)
# data <- rnorm(100)
# result <- nppi_optimal_block_length(data, stat_function)
# print(result)
