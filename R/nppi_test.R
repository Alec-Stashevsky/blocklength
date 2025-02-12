# Load necessary package
if (!requireNamespace("tseries", quietly = TRUE)) {
  stop("The 'tseries' package is required but not installed.")
}

#' Nonparametric Plug-in (NPPI) Block Length Selection
#' @param data Time series data
#' @param statistic Function to compute statistic of interest
#' @param initial_block_length Initial guess for block length (defaults to n^(1/3))
#' @param num_samples Number of bootstrap samples
#' @param d Exponent for bias correction, default is 1 (may depend on statistic)
#' @param ... Additional arguments for the statistic function
#' @return Optimal block length
nppi <- function(data, statistic, initial_block_length = NULL, num_samples = 1000, ...) {
  n <- length(data)

  # Step 1: Initial Block Length Selection
  if (is.null(initial_block_length)) {
    initial_block_length <- floor(n^(1/5))  # Matches theory in Lahiri et al. (2007)
  }

  # Step 2: Bias Estimation using block bootstrap
  l <- initial_block_length
  bootstrap_results_l <- tseries::tsbootstrap(data, nb = num_samples, b = l, type = "block")
  bootstrap_results_2l <- tseries::tsbootstrap(data, nb = num_samples, b = 2 * l, type = "block")

  # Here we assume the statistic is the mean,
  est_l <- mean(apply(bootstrap_results_l, 2, statistic, ...))
  est_2l <- mean(apply(bootstrap_results_2l, 2, statistic, ...))

  bias_estimate <- 2 * (est_l - est_2l)  # Corrected Bias formula

  # Step 3: Variance Estimation using JAB
  jab_variance <- jab_variance_estimate(data, statistic, l, num_samples, ...)

  # Step 4: Optimal Block Length Calculation (Formula-Based)
  optimal_block_length <- floor((2 * bias_estimate^2 / (1 * jab_variance))^(1 / 3) * n^(1 / 3))

  return(optimal_block_length)
}



#' Jackknife-After-Bootstrap (JAB) Variance Estimation
#' @param data Time series data
#' @param statistic Function to compute statistic of interest
#' @param block_length Block size for the bootstrap
#' @param num_samples Number of bootstrap samples
#' @param ... Additional arguments for the statistic function
#' @return Jackknife variance estimate
jab_variance_estimate <- function(data, statistic, block_length, num_samples = 1000, ...) {
  n <- length(data)

  # Optimal block deletion size
  m <- max(1, floor(n^(1/3) * block_length^(2/3)))

  jackknife_estimates <- numeric(n - m)

  # Compute the original bootstrap estimate
  bootstrap_results <- tseries::tsbootstrap(data, nb = num_samples, b = block_length, type = "block")
  original_estimate <- mean(apply(bootstrap_results, 2, statistic, ...))

  # Perform Jackknife-After-Bootstrap (JAB) with blocks of blocks deletion
  for (i in 1:(n - m)) {
    jackknife_sample <- data[-(i:(i + m - 1))]  # Remove a block of size m
    bootstrap_results <- tseries::tsbootstrap(jackknife_sample, nb = num_samples, b = block_length, type = "block")
    jackknife_estimates[i] <- mean(apply(bootstrap_results, 2, statistic, ...))
  }

  # Calculate Jackknife variance estimate
  jackknife_mean <- mean(jackknife_estimates)
  jab_variance <- ((n - m) / m) * mean((jackknife_estimates - jackknife_mean)^2)

  return(jab_variance)
}
