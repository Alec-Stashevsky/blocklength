
#' HHJ Algorithm
#'
#' Perform the Hall, Horowitz, and Jing (1995) "HHJ" algorithm to select the
#' optimal block-length \eqn{(l)} for a block bootstrap.
#'
#' @param series a numeric vector or time series giving the original data for
#'  which to find the optimal block length for.
#' @param nb number of bootstrap series to compute.
#' @param n.iter maximum number of iterations for HHJ algorithm.
#' @param pilot.block.length pilot block length (\eqn{l*} \emph{in HHJ})
#'  for which to perform initial block bootstraps.
#' @param sub.block.size length of each overlapping subsample
#'  (\eqn{m} \emph{in HHJ}).
#' @param bofb length of the basic blocks in the \emph{block of blocks}
#'  bootstrap.
#' @param search.grid the range of solutions around l* to evaluate within the
#'  MSE function after 1st iteration.
#' @param grid.step number to increment over subsample block lengths.
#'  If grid.step = 1 then each block length will be evaluated in the MSE
#'  function, if grid.step > 1, the the MSE function will search over the
#'  sequence of block lengths from 1 to m by grid.step. If grid.step is supplied
#'  as a vector of length 2, the the first iteration will step by the first
#'  element and subsequent iterations will step by the second element.
#' @param cl a cluster object, created by package \pkg{parallel} or by
#'  package \pkg{snow}. If \code{NULL}, use the non-parallel method.
#' @param verbose a logical value, if set to \code{FALSE} then no interim
#'  messages are output to the console. Error messages will still be output.
#'  Default is \code{TRUE}.
#' @param plots a logical value, if set to \code{FALSE} then no interim
#'  plots are output to the console. Default is \code{TRUE}.
#'
#' @export
#'
#' @examples
#' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' # Calculate optimal block length for series
#' hhjboot(sim, sub.block.size = 10)
#'
#' \dontrun{
#' # Use parallel computing
#' library(parallel)
#'
#' # Make cluster object with all cores available
#' cl <- makeCluster(detectCores())
#'
#' # Calculate optimal block length for series
#' hhjboot(sim, cl = cl)
#' }
#'
hhjboot <- function(series,
                    nb = 100L,
                    n.iter = 10L,
                    pilot.block.length = NULL,
                    sub.block.size = NULL,
                    bofb = 1L,
                    search.grid = NULL,
                    grid.step = c(1L, 1L),
                    cl = NULL,
                    verbose = TRUE,
                    plots = TRUE) {

  # Check arguments
  stopifnot(class(series) %in% c("integer", "numeric", "ts"))
  stopifnot(all.equal(nb %% 1, 0))
  stopifnot(all.equal(n.iter %% 1, 0))
  stopifnot(class(grid.step) %in% c("integer", "numeric"))

  if (length(grid.step) < 2) {
    stopifnot(all.equal(grid.step %% 1, 0))
    grid.step <- c(grid.step, grid.step)
  } else if (length(grid.step) > 2) {
    stop("grid.step must be a vector of at most length 2")
  }

  if (!is.null(search.grid)) {
    stopifnot(all.equal(search.grid %% 1, 0))
  }

  # Save function call
  call <- match.call()

  # Length of whole series
  n <- length(series)

  # Set pilot block-length
  if (is.null(pilot.block.length)) {
    l_star <- round(n^(1 / 5))
  } else {
    l_star <- round(pilot.block.length)
  }

  # Print pilot message
  if (isTRUE(verbose)) {
    message(" Pilot block length is: ", l_star)
  }

  # Check block-length of subsamples
  if (is.null(sub.block.size)) {
    m <- round(n^(1 / 5) * n^(1 / 3))
  } else {
    m <- sub.block.size
  }

  # Initialize overlapping sub samples list
  series.list <- vector(mode = "list", length = length((n - m + 1)))

  for (j in 1:n.iter) {

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

    # Search total grid on first iteration then +/- search.grid over next
    if (j == 1) {
      search_grid <- list(seq(from = 1, to = m, by = grid.step[1]))
    } else if (!is.null(search.grid)) {
      search_grid <- list(seq(
        from = max(1, l_m - search.grid),
        to = min(m, l_m + search.grid),
        by = grid.step[2]
      ))
    } else {
      search_grid <- list(seq(from = 1, to = m, by = grid.step[2]))
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
        lengths(search_grid), " function evaluations required.")
    }

    for (i in seq_len(length.out = (n - m + 1))) {

      # Create sub-blocks of length m = sub.block.size
      series.list[[i]] <- series[seq(from = i, to = (i + m - 1), by = 1)]
    }

    # Run optimization process
    if (is.null(cl)) {

      # Optimize MSE function over l in serial
      sol <- sapply(X = search_grid[[1]], FUN = hhjMSE,
        # Bootstrap parameters for hhjMSE
        n = n,
        m = m,
        series.list = series.list,
        nb = nb,
        bofb = bofb,
        v_star = v_star)

    } else {

      # Optimize MSE function over l in parallel
      sol <- parallel::parSapply(cl = cl, X = search_grid[[1]], FUN = hhjMSE,
        # Bootstrap parameters for hhjMSE
        n = n,
        m = m,
        series.list = series.list,
        nb = nb,
        bofb = bofb,
        v_star = v_star)
    }

    # Save plot data
    p.data <- data.frame(
      Iteration = j,
      Grid = round(search_grid[[1]] * ((n / m)^(1 / 3))),
      MSE = sol
      )

    # Plot MSE over l and color minimizing value red
    if (isTRUE(plots)) {
      plot(
        x = p.data$Grid,
        y = p.data$MSE,
        main = paste0(
          "MSE Plot for: ",
          deparse(substitute(series)),
          "\n",
          "Iteration: ",
          j
        ),

        xlab = "Block Length (l)",
        ylab = "MSE",
        col = ifelse(p.data$MSE == min(p.data$MSE), "red", "black")
      )
    }

    # Save l that minimizes MSE of subsample blocks
    l_m <- round(which.min(sol))

    # Break if l_m converges to previous l*
    if (l_star == round((n / m)^(1 / 3) * l_m)) {

      # Print convergence message
      if (isTRUE(verbose)) {
        message(" Converged at block length (l): ", round(l_star))
      }

      # Compile results list with custom class
      result <- structure(list(
        "Optimal Block Length" = l_star,
        "Subsample block size (m)" = m,
        "MSE Data" = p.data,
        "Call" = call
        ), class = "hhjboot")

      # Return list of results
      return(result)
    }

    # Scale l_m back to l_n for next iteration
    l_star <- round(((n / m)^(1 / 3)) * l_m)

    # Print iteration message
    if (isTRUE(verbose)) {
      message(" Chosen block length: ", l_star, "  After iteration: ", j)
    }
  }
}

# Helper Functions --------------------------------------------------------

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
