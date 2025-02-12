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




# quick test --------------------------------------------------------------

data <- c(-0.469391281673474, 0.343097144910637, 0.186715750040433, -0.531079700083545,
          -0.00678837122112312, 0.272434739497406, -0.252716426750953,
          -0.628788348685316, 3.47718937839464, 4.57651913366722, 0.296657442215465,
          -0.74850868705551, 2.14110624965374, 2.05112503816298, 1.67298490830627,
          1.90661318089383, -1.1727613558233, 0.238819327615459, 3.45874520333373,
          1.8405751291312, -1.30513889465471, -1.31563535093195, -1.26184834922138,
          -1.18694295176143, -1.28844555000195, -1.38758387043449, -0.85227121184701,
          0.715103155786832, 0.36789114896178, -1.1921918662436, -1.09197587115688,
          0.40812049159522, 0.0951237145441158, -0.365248847235551, -0.315721717061959,
          -1.26232128514309, -1.02079522254454, 1.41306990392656, 1.80590718861899,
          -0.644972335369019, -1.31519298934983, -0.798868334814169, -0.379002846958051,
          -0.797396341877188, 1.08377542127941, 1.69438404490305, -0.390370535540874,
          -1.15651782493652, -1.37617572695952, -1.03537947492713, -0.328494359114517,
          -0.693582633148482, -1.4085150217723, -1.21804600551757, -0.776457912635603,
          -0.929522879686841, -1.21868873908037, -0.908347984012241, -0.485022563679374,
          -0.663758166431915, -0.380673249388337, -0.422387998507367, -1.24650293141019,
          2.0479568491532, 2.04606311420025, -0.886927325252736, -0.919638488811108,
          -0.714895250046149, -0.689972850376533, -1.364284359842, -1.18727578710894,
          -1.20786682539434, -1.25769853823518, -0.804184090047034, -0.853551660898697,
          -1.27353402117699, -1.20784213199385, -1.19846482087295, -0.00388516673087015,
          1.1123105386246, -0.2220642752409, -0.670667690018065, -0.703222090336591,
          -1.33536129811158, -1.33419831139096, -1.08156833154194, -0.87693993904771,
          0.769398500868643, 1.26218438757083, 0.0857996813309342, -0.567148833936639,
          -0.0956653870660815, 3.20823245486653, 1.94350856334395, -0.782005104017045,
          -0.789494120012539, -1.22190777052684, -0.296135543929849, -0.209322831007107,
          -0.444213158260646, -0.613649061906346, -1.24112843793464, -0.429059054611025,
          0.352362935921961, 1.53964060250313, 1.4434480717644, -0.51388118032269,
          -0.86505886543298, -0.596854823362347, -0.992142717269907, -1.40952942586341,
          -1.08624700497268, 0.315698488228275, 0.920593870504548, -0.0603815350192414,
          -0.815576619925301, -0.948299811580347, -0.896783229056686, 0.433568547394551,
          0.809685706136017, -0.72332072411797, -1.02338145274977, -1.10187159087108,
          -1.01662571754982, -1.01840956407334)

result <- nppi(data, stat_function = mean, num_bootstrap = 300)
result


# Simulation Study  -------------------------------------------------------

# Define the model (6.1): X_t = (epsilon_t + epsilon_{t-1}) / sqrt(2), epsilon_t ~ chi^2(1) - 1
generate_data_model_6_1 <- function(n) {
  epsilon <- rchisq(n + 1, df = 1) - 1
  X <- (epsilon[-1] + epsilon[-(n + 1)]) / sqrt(2)
  return(X)
}

# Run simulation to estimate optimal block length
simulate_optimal_block_length <- function(num_simulations = 2, n = 125, block_lengths = 2:4) {
  mse_results <- numeric(length(block_lengths))


  for (l in block_lengths) {
    print(paste("Block Length:", l))
    pb <- txtProgressBar(min = 0, max = num_simulations, style = 3)
    estimates <- numeric(num_simulations)

    for (sim in 1:num_simulations) {
      data <- generate_data_model_6_1(n)
      result <- nppi(data, stat_function = function(x) mean(x), initial_block_size = l, num_bootstrap = 300)
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
