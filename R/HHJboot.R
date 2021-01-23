
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
#'
#' @export
#'
#' @examples
#'
#' # Generate AR(1) time series
#'sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'   n = 500, innov = rnorm(500))
#'
#' # Calculate optimal block length for series
#' hhjboot(sim)
#'
#' # Only evaluate every other block length after 1st iteration
#' hhjboot(sim, grid.step = 2)
#'
#' # Only evaluate +/- 20 block length from 1st iteration's solution
#' hhjboot(sim, search.grid = 20)
#'
#' \dontrun{
#'  # Use parallel computing
#'  library(parallel)
#'
#'  # Make cluster object with all cores available
#'  cl <- makeCluster(detectCores())
#'
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
                    cl = NULL) {

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

  # Function argument tests
  if (is.null(pilot.block.length)) {
    l_star <- round(n^(1 / 5))
  } else {
    l_star <- round(pilot.block.length)
  }

  # Print pilot message
  message(" Pilot block length is: ", l_star)

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


    if (!is.null(cl)) {

      # Send internal MSE function parameters to each cluster
      parallel::clusterExport(cl = cl, list(
        "n", "m", "series.list", "nb", "bofb",
        "v_star", "l_star"
      ), envir = environment())
    }


    if (j == 1) {

      # Output message to wait for gridSearch to process
      message("Performing minimization may take some time")
      message(
        "Calculating MSE for each level in subsample: ",
        lengths(search_grid), " function evaluations required."
      )
    }

    for (i in 1:(n - m + 1)) {

      # Create sub-blocks of length m = sub.block.size
      series.list[[i]] <- series[i:(i + m - 1)]
    }

    if (is.null(cl)) {

      # Optimize MSE function over l
      sol <- sapply(X = search_grid[[1]], FUN = function(l) {

        # Initialize Squared Error Vector
        se <- rep(NA, (n - m + 1))

        # Bootstrap each subsample, non-parallel computation
        output <- lapply(series.list, tseries::tsbootstrap,
          statistic = stats::var,
          type = "block",
          nb = nb,
          b = l,
          m = bofb
        )

        for (i in 1:length(output)) {

          # Calculate Squared Error
          se[i] <- mean(output[[i]]$statistic - v_star)^2
        }

        # Calculate MSE
        return(mean(se))
      })
    } else {

      # Optimize MSE function over l
      sol <- parallel::parSapply(cl = cl, X = search_grid[[1]], FUN = function(l) {


        # Initialize Squared Error and Bootstrap Vector
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

        for (i in 1:length(output)) {

          # Calculate Squared Error
          se[i] <- mean(output[[i]]$statistic - v_star)^2
        }

        # Calculate MSE
        return(mean(se))
      })
    }

    # Plot MSE over l and color minimizing value red
    plot(
      x = search_grid[[1]] * ((n / m)^(1 / 3)), y = sol,
      main = paste0(
        "MSE Plot for: ", deparse(substitute(series)), "\n",
        "Iteration: ", j
      ),
      xlab = "Block Length (l)", ylab = "MSE",
      col = ifelse(sol == min(sol), "red", "black")
    )

    # Save l that minimizes MSE of subsample blocks
    l_m <- round(which.min(sol))

    # Break if l_m equal to previous l*
    if (l_star == round((n / m)^(1 / 3) * l_m)) {

      # Print convergence message
      message(" Converged at block length (l): ", round(l_star))

      # Compile Result list
      result <- list(
        "Optimal Block Length" = l_star,
        "Subsample block size (m)" = m,
        "Call" = call
      )

      # Return list of results
      return(result)
    }

    # Scale l_m back to l_n for next iteration
    l_star <- round(((n / m)^(1 / 3)) * l_m)

    # Print iteration message
    message(" Chosen block length: ", l_star, "  After iteration: ", j)
  }
}
