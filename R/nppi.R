# NPPI Algorithm for Optimal Block Length in Block Bootstrap

# Function to estimate optimal block length using NPPI algorithm
nppi <- function(data, stat_function, r = 1, a = 1, initial_block_size = NULL, num_bootstrap = 1000, c_1 = 1L) {
  n <- length(data)

  # Step 1: Set initial block size if not provided
  if (is.null(initial_block_size)) {
    initial_block_size <- max(2, round(c_1*n^(1 / (r + 4))))
  }

  # Set c_2 to 1 if r=1 and to 0.1 otherwise
  c_2 <- ifelse(r == 1, 1, 0.1)

  # Function to create overlapping blocks
  create_overlapping_blocks <- function(data, block_size) {
    if (block_size > length(data)) {
      stop("Block size is larger than data length.")
    }

    N <- length(data) - block_size + 1
    blocks <- embed(data, block_size)[N:1, ]

    if (is.null(dim(blocks))) {
      blocks <- matrix(blocks, ncol = block_size)
    }

    return(blocks)
  }

  # Step 1: Manual Moving Block Bootstrap (MBB) Resampling
  mbb_resample <- function(blocks, num_samples) {
    N <- nrow(blocks)

    resampled_blocks <- replicate(num_samples, {
      sampled_indices <- sample(1:N, N, replace = TRUE)
      sampled_blocks <- blocks[sampled_indices, , drop = FALSE]  # Keep as matrix
      return(sampled_blocks)
    }, simplify = FALSE)

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
  C_hat_1 <- n * initial_block_size^(-r) * variance
  C_hat_2 <- initial_block_size * bias

  optimal_block_length <- ((2 * C_hat_2^2) / (r * C_hat_1)^(1 / (r + 2)) * n^(1 / (r + 2)))

  return(list(optimal_block_length = optimal_block_length, bias = bias, variance = variance))

}



# Simulation Study  -------------------------------------------------------

# Define the model (6.1): X_t = (epsilon_t + epsilon_{t-1}) / sqrt(2), epsilon_t ~ chi^2(1) - 1
generate_data_model_6_1 <- function(n) {
  epsilon <- rchisq(n + 1, df = 1) - 1
  X <- (epsilon[-1] + epsilon[-(n + 1)]) / sqrt(2)
  return(X)
}

# Run simulation to estimate optimal block length
simulate_optimal_block_length <- function(num_simulations = 500, n = 125, block_lengths = 2:10) {
  mse_results <- numeric(length(block_lengths))


  for (l in block_lengths) {
    print("Running simulation for block length: ", l)
    pb <- txtProgressBar(min = 0, max = num_simulations, style = 3)
    estimates <- numeric(num_simulations)

    for (sim in 1:num_simulations) {
      data <- generate_data_model_6_1(n)
      result <- nppi(data, stat_function = function(x) mean(x), initial_block_size = l)
      estimates[sim] <- result$optimal_block_length
      setTxtProgressBar(pb, sim)
    }

    mse_results[l] <- mean((estimates - 3)^2)  # Compare with true optimal block length = 3
  }

  return(data.frame(Block_Length = block_lengths, MSE = mse_results))
}

# Run the simulation
set.seed(32)
data <- generate_data_model_6_1(125)
simulation_results <- simulate_optimal_block_length()
print(simulation_results)
