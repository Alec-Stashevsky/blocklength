
#' Title
#'
#' @param series
#' @param statistic
#' @param nb
#' @param n.iter
#' @param pilot.block.length
#' @param subsample.size
#' @param bofb
#' @param search.grid
#'
#' @return
#' @export
#'
#' @examples
#'
#'
optim.HHJ.21 <- function(series, nb = 100L, n.iter = 10,
                pilot.block.length = NULL,
                sub.block.size = NULL,
                bofb = 1,
                search.grid = NULL,
                cl = NULL) {

  call <- match.call()

  # Save n of series
  n <- length(series)

  if (is.null(pilot.block.length)) {

    l_star <- round(n^(1/5))

  } else { l_star <- pilot.block.length }

  if (is.null(sub.block.size)) {

    m <- round(n^(1/5)*n^(1/3))

  } else { m <- sub.block.size }

  if (is.null(cl)){

    grid.method <- NULL

  } else { grid.method <- 'snow' }

  # Print pilot message
  message(" Pilot block length is: ", l_star)

  for (j in 1:n.iter){

    # Bootstrap whole series
    boot_temp <- tseries::tsbootstrap(series,
                                      statistic = var,
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

    if (j == 1) {

    # Output message to wait for gridSearch to process
    message("Performing minimization may take some time")

      }

    # Optimize MSE function over l
    sol <- NMOF::gridSearch(function(l) {

      # Initialize Squared Error and Bootstrap Vector
      se <- rep(NA, (n - m + 1))
      series.list <- vector(mode = "list", length = length(se))

      for (i in 1:(n - m + 1)) {
        # Create sub-blocks of length m = sub.block.size
        series.list[[i]] <- series[i:(i + m - 1)]

      }


      # if (is.null(cl)){

        # Use non-parallel computation
        output <- lapply(series.list, tseries::tsbootstrap,
                         statistic = var,
                         type = "block",
                         nb = nb,
                         b = l,
                         m = bofb)

      # } else {
      #
      #   # Use load balanced parallelization
      #   output <- parallel::parLapplyLB(cl = cl, X = series.list,
      #                                   fun = tseries::tsbootstrap,
      #                                   statistic = var,
      #                                   type = "block",
      #                                   nb = nb,
      #                                   b = l,
      #                                   m = bofb)
      # }


      for (i in 1:length(output)){

        se[i] <- mean(output[[i]]$statistic - v_star)^2

      }

      # Calculate MSE
      return(mean(se))

      },

      levels = search_grid,
      method = 'snow',
      cl = cl)

    # Plot MSE over l
    plot(x = (1:m)*((n/m)^(1/3)), y = sol$values,
         main = paste0("Iteration: ", j),
         xlab = "Block Length (l)", ylab = "MSE")

    # Save l that minimizes MSE
    l_m <- round(sol$minlevels)

    # Break if l_m equal to previous l*
    if (l_star == round((n/m)^(1/3)*l_m)) {

      # Stop cluster
      parallel::stopCluster(cl = cl)

      # Print convergence message
      message(" Converged at block length (l): ", round(l_star))

      # Compile Result list
      result <- list("Optimal Block Length" = l_star,
                     "Pilot Number of Blocks (m)" = m,
                     "Call" = call)

      return(result)

    }

    # Scale l_m back to l_n for next iteration
    l_star <- round((n/m)^(1/3)*l_m)

    # Print iteration message
    message(" Chosen block length: ", l_star, "  After iteration: ", j)

  }

}


cl <- makePSOCKcluster(3)
system.time(optim.HHJ.21(sim, cl = cl))
cl <- makePSOCKcluster(3)
system.time(optim.HHJ.2(sim, cl = cl))
