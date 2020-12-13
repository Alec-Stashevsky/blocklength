
#' HHJ Algorithm
#'
#' Perform the HHJ algorithm to select the optimal block length (l) for the
#' moving block bootstrap.
#'
#' @param series
#' @param nb
#' @param n.iter
#' @param pilot.block.length
#' @param sub.block.size
#' @param bofb
#' @param search.grid
#' @param cl
#'
#' @return
#' @export
#'
#' @examples
#'
#' HHJboot(arima.sim(list(order = c(1, 0, 0), ar = 0.5, n = 500, innov = rnorm(500))))
#'
HHJboot <- function(series, nb = 100L, n.iter = 10,
                    pilot.block.length = NULL,
                    sub.block.size = NULL,
                    bofb = 1,
                    search.grid = NULL,
                    cl = NULL) {

  # Save function call
  call <- match.call()

  # Length of whole series
  n <- length(series)

  # Function argument tests
  if (is.null(pilot.block.length)) {

    l_star <- round(n^(1/5))

  } else { l_star <- pilot.block.length }

  # Print pilot message
  message(" Pilot block length is: ", l_star)

  if (is.null(sub.block.size)) {

    m <- round(n^(1/5)*n^(1/3))

  } else { m <- sub.block.size }

  # Initialize overlapping sub samples list
  series.list <- vector(mode = "list", length = length((n - m + 1)))

  for (j in 1:n.iter){

    # Bootstrap variance of whole series
    boot_temp <- tseries::tsbootstrap(series,
                                      statistic = stats::var,
                                      type = "block",
                                      nb = nb,
                                      b = l_star,
                                      m = bofb)

    # Save updated variance of whole series
    v_star <- mean(boot_temp$statistic)

    # Search total grid on first iteration then +/- search.grid over next
    if (j == 1) {

      search_grid <- list(1:m)

      } else if (is.integer(search.grid)) {

        search_grid <- list(max(1, l_m - search.grid):min(m, l_m + search.grid))

        } else {

          search_grid <- list(1:m)

        }


    if (!is.null(cl)){

      # Send internal MSE function parameters to each cluster
      clusterExport(cl = cl, list("n", "m", "series.list", "nb", "bofb",
                                "v_star", "l_star"), envir = environment())
    }


    if (j == 1) {

      # Output message to wait for gridSearch to process
      message("Performing minimization may take some time")
      message("Calculating MSE for each level in subsample: ",
              lengths(search_grid), " function evaluations required.")

    }

    for (i in 1:(n - m + 1)) {

      # Create sub-blocks of length m = sub.block.size
      series.list[[i]] <- series[i:(i + m - 1)]

    }

    if (is.null(cl)){

      # Optimize MSE function over l
      sol <- vapply(X = 1:m, FUN = function(l) {

        # Initialize Squared Error Vector
        se <- rep(NA, (n - m + 1))

        # Bootstrap each subsample, non-parallel computation
        output <- lapply(series.list, tseries::tsbootstrap,
                         statistic = stats::var,
                         type = "block",
                         nb = nb,
                         b = l,
                         m = bofb)

        for (i in 1:length(output)){

          # Calculate Squared Error
          se[i] <- mean(output[[i]]$statistic - v_star)^2

        }

        # Calculate MSE
        return(mean(se))

      })

    } else {

      # Optimize MSE function over l
      sol <- parallel::parSapply(cl = cl, X = 1:m, FUN = function(l) {


         # Initialize Squared Error and Bootstrap Vector
         se <- rep(NA, (n - m + 1))

         # Bootstrap each subsample, non-parallel computation
         output <- lapply(X = series.list,
                          FUN = tseries::tsbootstrap,
                          statistic = var,
                          type = "block",
                          nb = nb,
                          b = l,
                          m = bofb)

         for (i in 1:length(output)){

           # Calculate Squared Error
           se[i] <- mean(output[[i]]$statistic - v_star)^2

         }

         # Calculate MSE
         return(mean(se))

       })
    }

    # Plot MSE over l and color minimizing value red
    plot(x = (1:m)*((n/m)^(1/3)), y = sol,
         main = paste0("MSE Plot for: ", deparse(substitute(series)), "\n",
                       "Iteration: ", j),
         xlab = "Block Length (l)", ylab = "MSE",
         col = ifelse(sol == min(sol), "red", "black"))

    # Save l that minimizes MSE of subsample blocks
    l_m <- round(which.min(sol))

    # Break if l_m equal to previous l*
    if (l_star == round((n/m)^(1/3)*l_m)) {

      # Print convergence message
      message(" Converged at block length (l): ", round(l_star))

      # Compile Result list
      result <- list("Optimal Block Length" = l_star,
                     "Pilot Number of Blocks (m)" = m,
                     "Call" = call)

      # Return list of results
      return(result)

    }

    # Scale l_m back to l_n for next iteration
    l_star <- round(((n/m)^(1/3))*l_m)

    # Print iteration message
    message(" Chosen block length: ", l_star, "  After iteration: ", j)

  }

}
