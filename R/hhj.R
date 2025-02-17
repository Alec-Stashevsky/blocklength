#' Hall, Horowitz, and Jing (1995) "HHJ" Algorithm to Select the Optimal Block-Length
#'
#' Perform the Hall, Horowitz, and Jing (1995) "HHJ" cross-validation algorithm
#' to select the optimal block-length for a bootstrap on dependent data
#' (block-bootstrap). Dependent data such as stationary time series are suitable
#' for usage with the HHJ algorithm.
#'
#' The HHJ algorithm is computationally intensive as it relies on a
#' cross-validation process using a type of subsampling to estimate the mean
#' squared error (\eqn{MSE}) incurred by the bootstrap at various block-lengths.
#'
#' Under-the-hood, \code{hhj()} makes use of \code{\link[tseries]{tsbootstrap}},
#' \emph{see} Trapletti and Hornik (2020), to perform the moving block-bootstrap
#' (or the \emph{block-of-blocks} bootstrap by setting \code{bofb > 1}) according
#' to Kunsch (1989).
#'
#' @param series a numeric vector or time series giving the original data for
#'  which to find the optimal block-length for.
#' @param nb an integer value, number of bootstrapped series to compute.
#' @param n_iter an integer value, maximum number of iterations for the HHJ
#'  algorithm to compute.
#' @param pilot_block_length a numeric value, the block-length (\eqn{l*} in
#'  \emph{HHJ}) for which to perform initial block bootstraps.
#' @param sub_sample a numeric value, the length of each overlapping
#'  subsample, \eqn{m} in \emph{HHJ}.
#' @param k a character string, either \code{"bias/variance"},
#'  \code{"one-sided"}, or \code{"two-sided"} depending on the desired object of
#'  estimation. If the desired bootstrap statistic is bias or variance then
#'  select \code{"bias/variance"} which sets \eqn{k = 3} per HHJ. If the object
#'  of estimation is the one-sided or two-sided distribution function, then set
#'  \code{k = "one-sided"} or \code{k = "two-sided"} which sets \eqn{k = 4} and
#'  \eqn{k = 5}, respectively. For the purpose of generating symmetric confidence
#'  intervals around an unknown parameter, \code{k = "two-sided"} (the default)
#'  should be used.
#' @param bofb a numeric value, length of the basic blocks in the
#'  \emph{block-of-blocks} bootstrap, \emph{see} \code{m = } for
#'  \code{\link[tseries]{tsbootstrap}} and Kunsch (1989).
#' @param search_grid a numeric value, the range of solutions around \eqn{l*} to
#'  evaluate within the \eqn{MSE} function \emph{after} the first iteration. The
#'  first iteration will search through all the possible block-lengths unless
#'  specified in \code{grid_step = }.
#' @param grid_step a numeric value or vector of at most length 2, the number of
#'  steps to increment over the subsample block-lengths when evaluating the
#'  \eqn{MSE} function. If \code{grid_step = 1} then each block-length will be
#'  evaluated in the \eqn{MSE} function. If \code{grid_step > 1}, the \eqn{MSE}
#'  function will search over the sequence of block-lengths from \code{1} to
#'  \code{m} by \code{grid_step}. If \code{grid_step} is a vector of length 2,
#'  the first iteration will step by the first element of \code{grid_step} and
#'  subsequent iterations will step by the second element.
#' @param cl a cluster object, created by package \pkg{parallel},
#'  \pkg{doParallel}, or \pkg{snow}. If \code{NULL}, no parallelization will be
#'  used.
#' @param verbose a logical value, if set to \code{FALSE} then no interim
#'  messages are output to the console. Error messages will still be output.
#'  Default is \code{TRUE}.
#' @param plots a logical value, if set to \code{FALSE} then no interim plots
#'  are output to the console. Default is \code{TRUE}.
#'
#' @return an object of class 'hhj'
#'
#' @section References:
#'
#' Adrian Trapletti and Kurt Hornik (2020). tseries: Time Series Analysis and
#'      Computational Finance. R package version 0.10-48.
#'
#' Kunsch, H. (1989) The Jackknife and the Bootstrap for General Stationary
#'      Observations. The Annals of Statistics, 17(3), 1217-1241. Retrieved
#'      February 16, 2021, from \doi{10.1214/aos/1176347265}
#'
#' Peter Hall, Joel L. Horowitz, Bing-Yi Jing, On blocking rules for the
#'      bootstrap with dependent data, Biometrika, Volume 82, Issue 3,
#'      September 1995, Pages 561-574, DOI: \doi{10.1093/biomet/82.3.561}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' # Calculate optimal block length for series
#' hhj(sim, sub_sample = 10)
#'
#'
#' # Use parallel computing
#' library(parallel)
#'
#' # Make cluster object with 2 cores
#' cl <- makeCluster(2)
#'
#' # Calculate optimal block length for series
#' hhj(sim, cl = cl)
#' }
#'
hhj <- function(
  series,
  nb = 100L,
  n_iter = 10L,
  pilot_block_length = NULL,
  sub_sample = NULL,
  k = "two-sided",
  bofb = 1L,
  search_grid = NULL,
  grid_step = c(1L, 1L),
  cl = NULL,
  verbose = TRUE,
  plots = TRUE
  ) {

  # Check arguments
  stopifnot(class(series) %in% c("integer", "numeric", "ts"))
  stopifnot(all.equal(nb %% 1, 0))
  stopifnot(all.equal(n_iter %% 1, 0))
  stopifnot(class(grid_step) %in% c("integer", "numeric"))

  if (length(grid_step) < 2) {
    stopifnot(all.equal(grid_step %% 1, 0))
    grid_step <- c(grid_step, grid_step)
  } else if (length(grid_step) > 2) {
    stop("grid_step must be a vector of at most length 2")
  }

  if (!is.null(search_grid)) {
    stopifnot(all.equal(search_grid %% 1, 0))
  }

  # Save function call
  call <- match.call()

  # Length of whole series
  n <- length(series)

  # Set pilot block-length
  if (is.null(pilot_block_length)) {
    l_star <- round(n^(1 / 5))
  } else {
    l_star <- round(pilot_block_length)
  }

  # Print pilot message
  if (isTRUE(verbose)) {
    message(" Pilot block length is: ", l_star)
  }

  # Set estimation type (k)
  if (k == "bias/variance") {
    k <- 3
  } else if (k == "one-sided") {
    k <- 4
  } else if (k == "two-sided") {
    k <- 5
  } else {
    stop("k must be in c('bias/variance', 'one-sided', 'two-sided')")
  }

  # Check block-length of subsamples
  if (is.null(sub_sample)) {
    m <- round(n^(1 / 5) * n^(1 / k))
  } else if (sub_sample >= n) {
    stop("sub_sample must be less than series length")
  } else {
    m <- round(sub_sample)
  }

  # Initialize overlapping sub samples list
  series.list <- vector(mode = "list", length = length((n - m + 1)))

  for (j in 1:n_iter) {

    # Bootstrap variance of whole series
    boot_temp <- tseries::tsbootstrap(series,
      statistic = stats::var,
      type = "block",
      nb = nb,
      b = l_star,
      m = bofb
    )

    # Save updated variance of whole series
    v_star <- mean(boot_temp$statistic)

    # Search total grid on first iteration then +/- search_grid over next
    if (j == 1) {
      grid <- list(seq(from = 1, to = m, by = grid_step[1]))
    } else if (!is.null(search_grid)) {
      grid <- list(seq(
        from = max(1, l_m - search_grid),
        to = min(m, l_m + search_grid),
        by = grid_step[2]
      ))
    } else {
      grid <- list(seq(from = 1, to = m, by = grid_step[2]))
    }

    # Setup environment for cluster workers
    if (!is.null(cl)) {

      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Package \"parallel\" needed for this function to work. Please install it.",
          call. = FALSE)
      }

      # Send internal MSE function parameters to each cluster
      parallel::clusterExport(cl = cl, list(
        "n", "m", "series.list", "nb", "bofb", "v_star"
        ), envir = environment())
    }

    # Minimization message
    if (j == 1 & isTRUE(verbose)) {

      # Output message to wait for gridSearch to process
      message("Performing minimization may take some time")
      message("Calculating MSE for each level in subsample: ",
        lengths(grid), " function evaluations required.")
    }

    # Create sub-blocks of length m = sub_sample
    for (i in seq_len(length.out = (n - m + 1))) {
      series.list[[i]] <- series[seq(from = i, to = (i + m - 1), by = 1)]
    }

    # Run optimization process
    if (is.null(cl)) {

      # Optimize MSE function over l in serial
      sol <- sapply(X = grid[[1]], FUN = hhjMSE,
        # Bootstrap parameters for hhjMSE
        n = n,
        m = m,
        series.list = series.list,
        nb = nb,
        bofb = bofb,
        v_star = v_star)

    } else {

      # Optimize MSE function over l in parallel
      sol <- parallel::parSapply(cl = cl, X = grid[[1]], FUN = hhjMSE,
        # Bootstrap parameters for hhjMSE
        n = n,
        m = m,
        series.list = series.list,
        nb = nb,
        bofb = bofb,
        v_star = v_star)
    }

    # Save plot data
    if (j == 1) {
      p.data <- data.frame(
        Iteration = j,
        BlockLength = round(grid[[1]] * ((n / m)^(1 / k))),
        MSE = sol
        )
    } else {
      p.data <- rbind(p.data, data.frame(
        Iteration = j,
        BlockLength = round(grid[[1]] * ((n / m)^(1 / k))),
        MSE = sol
        ))
    }

    # Plot MSE over l and color minimizing value red
    if (isTRUE(plots)) {
      plot(
        x = round(grid[[1]] * ((n / m)^(1 / k))),
        y = sol,
        main = paste0(
          "MSE Plot for: ",
          deparse(substitute(series)),
          "\n",
          "Iteration: ",
          j
        ),

        xlab = "Block Length (l)",
        ylab = "MSE",
        col = ifelse(sol == min(sol), "red", "black")
      )
    }

    # Save l that minimizes MSE of subsample blocks
    l_m <- round(which.min(sol))

    # Break if l_m converges to previous l*
    if (l_star == round((n / m)^(1 / k) * l_m)) {

      # Print convergence message
      if (isTRUE(verbose)) {
        message(" Converged at block length (l): ", round(l_star))
      }

      # Compile results list with custom class
      result <- structure(
        list(
          "Optimal Block Length" = l_star,
          "Subsample block size (m)" = m,
          "MSE Data" = p.data,
          "Iterations" = j,
          "Series" = deparse(substitute(series)),
          "Call" = call
          ),
        class = "hhj")

      # Return list of results
      return(result)
    }

    # Scale l_m back to l_n for next iteration
    l_star <- round(((n / m)^(1 / k)) * l_m)

    # Print iteration message
    if (isTRUE(verbose)) {
      message(" Chosen block length: ", l_star, "  After iteration: ", j)
    }

    # Print if final-iteration is reached w/ no convergence
    if (j == n_iter) {

      # Warning message about no convergence
      warning(
        "Block-length has not converged. Stopping at iteration limit: ", j, "\n",
        "  Result may be spurious"
        )

      # Compile results list with custom class
      result <- structure(
        list(
          "Optimal Block Length" = l_star,
          "Subsample block size (m)" = m,
          "MSE Data" = p.data,
          "Iterations" = j,
          "Series" = deparse(substitute(series)),
          "Call" = call
        ),
        class = "hhj")

      # Return list of results
      return(result)
    }

  }
}

# Helper Functions --------------------------------------------------------

# MSE function to optimize
hhjMSE <- function(l, n, m, series.list, nb, bofb, v_star) {

  # Initialize Squared Error Vector
  se <- rep(NA, (n - m + 1))

  # Bootstrap each subsample, non-parallel computation
  output <- lapply(
    X = series.list,
    FUN = tseries::tsbootstrap,
    statistic = stats::var,
    type = "block",
    nb = nb,
    b = l,
    m = bofb
  )

  for (i in seq_len(length.out = length(output))) {

    # Calculate Squared Error
    se[i] <- (mean(output[[i]]$statistic) - v_star)^2
  }

  # Calculate MSE
  return(mean(se))

}
