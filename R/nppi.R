#' Lahiri, Furukawa, Lee, (2007) Nonparametric Plug-In "NPPI" Rule to Select the Optimal Block-Length
#'
#' This function implements the Nonparametric Plug-In (NPPI) algorithm,
#'  as proposed by Lahiri, Furukawa, and Lee (2007), to select the optimal block
#'  length for block bootstrap procedures.
#'  The NPPI method estimates the optimal block length by balancing bias and
#'  variance in block bootstrap estimators, particularly for time series and
#'  other dependent data structures. The function also leverages the Moving
#'  Block Bootstrap (MBB) method of (Künsch, 1989) and the Moving Blocks
#'  Jackknifte (MBJ) of Liu and Singh (1992).
#'
#'  Jackknife-After-Bootstrap (JAB) variance estimation (Lahiri, 2002).
#'
#' @param data A numeric vector representing the time series or dependent data.
#' @param stat_function A function to compute the statistic of interest
#'  (*e.g.*, mean, variance). The function should accept a numeric vector as input
#'  and return a scalar value (default is \code{mean}).
#' @param r The rate parameter for the MSE expansion (default is 1). This
#'  parameter controls the convergence rate in the bias-variance trade-off.
#' @param a The bias exponent (default is 1). Adjust this based on the
#'  theoretical properties of the statistic being bootstrapped.
#' @param initial_block_size Optional. The initial block size for bias estimation.
#'   If not provided, it is set to \code{max(2, round(c_1 * n^{1 / (r + 4)}))}, where \code{n} is the sample size.
#' @param m Optional. The number of blocks to delete in the
#'  Jackknife-After-Bootstrap (JAB) variance estimation. If not provided,
#'  it defaults to \code{floor(c_2 * n^{1/3} * initial_block_size^{2/3})}.
#' @param num_bootstrap The number of bootstrap replications for bias estimation
#'  (default is 1000).
#' @param c_1 A tuning constant for initial block size calculation (default is 1).
#' @param epsilon A small constant added to the variance to prevent division by
#'  zero (default is \code{1e-8}).
#'
#' @return A object of class \code{nppi} with the following components:
#' \describe{
#'   \item{optimal_block_length}{The estimated optimal block length for the block bootstrap procedure.}
#'   \item{bias}{The estimated bias of the block bootstrap estimator.}
#'   \item{variance}{The estimated variance of the block bootstrap estimator using the JAB method.}
#' }
#'
#'  @section References:
#'
#'  Efron, B. (1992), 'Jackknife-after-bootstrap standard errors and influence
#'      functions (with discussion)', Journal of Royal Statistical Society,
#'      Series B 54, 83-111.
#'
#' Kunsch, H. (1989) The Jackknife and the Bootstrap for General Stationary
#'      Observations. The Annals of Statistics, 17(3), 1217-1241. Retrieved
#'      February 16, 2021, from \url{http://www.jstor.org/stable/2241719}
#'
#' Lahiri, S. N., Furukawa, K., & Lee, Y.-D. (2007). A nonparametric plug-in
#'      rule for selecting optimal block lengths for Block Bootstrap Methods.
#'      Statistical Methodology, 4(3), 292–321. DOI:
#'      \doi{10.1016/j.stamet.2006.08.002}
#'
#' Lahiri, S. N. (2003). 7.4 A Nonparametric Plug-in Method. In Resampling
#'      methods for dependent data (pp. 186–197). Springer.
#'
#' Liu, R. Y. and Singh, K. (1992), Moving blocks jackknife and bootstrap
#'      capture weak dependence, in R. Lepage and L. Billard, eds, 'Exploring
#'      the Limits of the Bootstrap', Wiley, New York, pp. 225-248.
#'

#' @export
#'
#' @examples
#' #' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' # Estimate the optimal block length for the sample mean
#' result <- nppi(data = sim, stat_function = mean, num_bootstrap = 500)
#'
#' print(result$optimal_block_length)
#'
nppi <- function(
    data,
    stat_function = mean,
    r = 1,
    a = 1,
    initial_block_size = NULL,
    m = NULL,
    num_bootstrap = 1000,
    c_1 = 1L,
    epsilon = 1e-8) {

  # Set parameters
  n <- length(data)

  if (is.null(initial_block_size)) {
    initial_block_size <- max(2, round(c_1 * n^(1 / (r + 4))))
  }

  # Recommended values
  c_2 <- if (r == 1) 1 else 0.1

  if (is.null(m)) {
    m <- floor(c_2 * n^(1/3) * initial_block_size^(2/3))
  }

  # Step 1: Estimate bias via block bootstrap
  bias_result <- bias_estimator(data, initial_block_size, stat_function, num_bootstrap)

  # Step 2: Estimate variance via moving blocks jackknife
  variance <- jab_variance_estimator(
    bias_result$resampled_blocks,
    bias_result$original_blocks,
    stat_function,
    m,
    bias_result$phi_hat_n
    )

  # Step 3: Compute the NPPI estimate for the optimal block length
  C_hat_1 <- n * initial_block_size^(-r) * variance
  C_hat_2 <- initial_block_size * bias_result$bias
  optimal_block_length <- (2 * C_hat_2^2 / ((r * C_hat_1)^(1 / (r + 2)))) * n^(1 / (r + 2))

  # Compile results with custom class
  result <- structure(
    list(
      optimal_block_length = optimal_block_length,
      bias = bias_result$bias,
      variance = variance
    ),
    class = "nppi"
  )

  return(result)

}


# Helper Functions --------------------------------------------------------

# Create overlapping blocks per MBB Method of Kunsch (1989)
create_overlapping_blocks <- function(data, block_size) {
  if (block_size > length(data)) {
    stop("Block size is larger than data length.")
  }
  N <- length(data) - block_size + 1

  # embed() creates blocks in reverse order so we re-order rows
  blocks <- stats::embed(data, block_size)[N:1, ]
  if (is.null(dim(blocks))) {
    blocks <- matrix(blocks, ncol = block_size)
  }
  return(blocks)
}


# Moving Block Bootstrap (MBB) Resampling
mbb <- function(blocks, num_samples) {
  N <- nrow(blocks)
  resampled_blocks <- replicate(num_samples, {
    sampled_indices <- sample(1:N, N, replace = TRUE)
    blocks[sampled_indices, , drop = FALSE]
  }, simplify = FALSE)

  return(resampled_blocks)
}


# Moving Block Bootstrap Bias Estimator
bias_estimator <- function(data, block_size, stat_function, num_bootstrap) {
  # Create blocks for block sizes l and 2l
  blocks_l1   <- create_overlapping_blocks(data, block_size)
  blocks_2l1  <- create_overlapping_blocks(data, 2 * block_size)

  # Generate bootstrap samples for each block set
  resampled_l1   <- mbb(blocks_l1, num_bootstrap)
  resampled_2l1  <- mbb(blocks_2l1, num_bootstrap)

  # Compute bootstrap replicates of the statistic
  theta_l1   <- sapply(resampled_l1, function(bs) stat_function(as.vector(t(bs))))
  theta_2l1  <- sapply(resampled_2l1, function(bs) stat_function(as.vector(t(bs))))

  # Compute the bias estimate: 2*(phi_hat(l1) - phi_hat(2l1))
  phi_hat_n_of_l <- mean(theta_l1)
  phi_hat_n_of_2l <- mean(theta_2l1)
  bias <- 2 * (phi_hat_n_of_l - phi_hat_n_of_2l)

  return(list(
    bias = bias,
    resampled_blocks = resampled_l1,
    original_blocks  = blocks_l1,
    phi_hat_n = phi_hat_n_of_l)
  )
}


# JAB Variance Estimator
jab_variance_estimator <- function(resampled_blocks, original_blocks, stat_function, m, phi_hat_n, epsilon = 1e-8) {
  N <- nrow(original_blocks)
  M <- N - m + 1
  K <- length(resampled_blocks)

  # Pre-compute row "signatures" for each bootstrap sample
  # Each bootstrap sample is a matrix; we collapse each row into a string
  # We do this to speed up the overlap check in the next step
  bs_signatures <- lapply(resampled_blocks, function(bs) {
    apply(bs, 1, paste, collapse = ",")
  })

  # Initialize vector to hold JAB point values
  jab_point_values <- numeric(M)

  # For each deletion block (of m consecutive original blocks)
  for (i in seq_len(M)) {
    # Select the m consecutive blocks to be "deleted"
    deletion_blocks <- original_blocks[i:(i + m - 1), , drop = FALSE]
    deletion_sigs <- apply(deletion_blocks, 1, paste, collapse = ",")

    # Identify bootstrap samples with no overlap
    non_overlap_indices <- which(sapply(bs_signatures, function(bs_sig) {
      !any(deletion_sigs %in% bs_sig)
    }))

    if (length(non_overlap_indices) > 0) {
      # Combine the non-overlapping bootstrap samples and compute the statistic
      combined_bs <- do.call(rbind, resampled_blocks[non_overlap_indices])
      jab_point_values[i] <- stat_function(as.vector(t(combined_bs)))
    } else {
      jab_point_values[i] <- NA  # No valid bootstrap samples for this deletion
    }
  }

  # Remove any NA values and compute the JAB variance
  jab_point_values <- jab_point_values[!is.na(jab_point_values)]
  jab_var <- (m / (N - m)) * (1 / M) * sum((jab_point_values - phi_hat_n)^2) + epsilon

  return(jab_var)
}
