# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

### HHJ Package

optim.HHJ <- function(series, statistic, nb = 100L, n.iter = 10,
                pilot.block.length = NULL,
                pilot.block.number = NULL,
                bofb = 1,
                search.grid = NULL) {

  call <- match.call()

  # Save n of series
  n <- length(series)

  if (is.null(pilot.block.length)) {

    l_star <- round(n^(1/5))

  } else { l_star <- pilot.block.length }

  if (is.null(pilot.block.number)) {

    m <- round(n^(1/5)*n^(1/3))

  } else { m <- pilot.block.number }

  # Print pilot message
  message(" Pilot block length is: ", l_star)

  for (j in 1:n.iter){

    # Bootstrap whole series
    boot_temp <- tseries::tsbootstrap(series,
                                      statistic = statistic,
                                      type = "block",
                                      nb = nb,
                                      b = l_star,
                                      m = bofb)

    # Save updated variance of n statistic
    v_star <- mean(boot_temp$statistic)

    # Search over total grid on first iteration, around +/- 20 over subsequent
    if (j == 1) {

      search_grid <- list(1:m)

    } else if (is.integer(search.grid)) {

      search_grid <- list(max(1, l_m - search.grid):min(m, l_m + search.grid))

    } else {

      search_grid <- list(1:m)

    }

    # Optimize MSE function over l
    sol <- NMOF::gridSearch(function(l) {

      # Initialize Squared Error Vector
      se <- rep(NA, (n - m + 1))

      for (i in 1:(n - m + 1)) {

        # Bootstrap subset blocks w/ length l
        boot_temp <- tseries::tsbootstrap(series[i:(i + m - 1)],
                                          statistic = statistic,
                                          type = "block",
                                          nb = nb,
                                          b = l,
                                          m = bofb)

        # Calculate Squared Error
        se[i] <- (mean(boot_temp$statistic) - v_star)^2

      }

      # Calculate MSE
      return(mean(se))

    },

    levels = search_grid)

    # Plot MSE over l
    plot(x = (1:m)*((n/m)^(1/3)), y = sol$values, main = paste0("Iteration: ", j),
         xlab = "Block Length (l)", ylab = "MSE")

    # Save l that minimizes MSE
    l_m <- round(sol$minlevels)

    # Break if l_m equal to previous l*
    if (l_star == round((n/m)^(1/3)*l_m)) {

      # Print convergence message
      message(" Converged at block length (l): ", round(l_star))

      result <- list("Optimal Block Length" = l_star,
                     "Pilot Number of Blocks (m)" = m,
                     "Call" = call)

      return(result)

    }

    # Update l* for next iteration
    l_star <- round((n/m)^(1/3)*l_m)

    # Print iteration message
    message(" Chosen block length: ", l_star, "  After iteration ", j)

  }

}
