#' Plot Stability of JAB Point Estimates
#'
#' S3 Method for objects of class 'nppi'
#' This function visualizes the JAB point estimates across deletion blocks
#'  indices used to estimate variance of the NPPI algorithm.
#'
#' @param x An object of class \code{nppi}, containing JAB point values.
#'
#' @inheritParams nppi
#' @inheritDotParams base::plot
#'
#' @return No return value, called for side effects
#'
#' @export
#'
#' @examples
#' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' # Estimate the optimal block length for the sample mean
#' result <- nppi(data = sim, stat_function = mean, num_bootstrap = 500, m =2)
#'
#' # Use s3 method
#' plot(result)
#'
plot.nppi <- function(x, ...) {
  if (is.null(x$jab_point_values) || all(is.na(x$jab_point_values))) {
    stop("Error: No valid JAB variance estimates available for plotting.")
  }

  # Define x-axis (deletion block indices)
  deletion_blocks <- seq_along(x$jab_point_values)

  # Compute mean variance for reference
  mean_variance <- mean(x$jab_point_values, na.rm = TRUE)

  # Plot variance estimates (dots and connecting line)
  plot(deletion_blocks, x$jab_point_values, type = "o", pch = 16, col = "black",
       xlab = "Deletion Block Index", ylab = "JAB Point Value",
       main = "Stability of JAB Point Values", ...)

  # Add horizontal dashed line for mean variance
  graphics::abline(h = mean_variance, col = "black", lty = 2)

  # Place legend outside the plot area (below the x-axis, centered)
  graphics::legend("topleft", horiz = FALSE,
         legend = c("JAB Point Value", "Mean"),
         col = "black", pch = c(16, NA), lty = c(NA, 2), lwd = 1, cex = 0.9)
}
