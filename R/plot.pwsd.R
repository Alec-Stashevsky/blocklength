#' Plot Correlogram for Politis and White Auto\eqn{-}Correlation Implied Hypothesis Test
#'
#' S3 Method for objects of class 'pwsd'
#' \emph{See} \code{?plot.acf} of the \pkg{stats} package for more customization
#'  options on the correlogram, from which \code{plot.pwsd} is based
#'
#' @param x an of object of class 'pwsd' or 'acf'
#' @param main an overall title for the plot, if no string is supplied a default
#'  title  will be populated. \emph{See} \code{\link[graphics]{title}}
#' @param ylim a numeric of length 2 giving the y-axis limits for the plot
#'
#' @inheritParams pwsd
#' @inheritDotParams base::plot
#'
#'
#' @return No return value, called for side effects
#'
#' @export
#'
#' @examples
#' # Use S3 Method
#'
#' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' b <- pwsd(sim, round = TRUE, correlogram = FALSE)
#' plot(b)
#'
plot.pwsd <- function(x, c = NULL, main = NULL, ylim = NULL, ...) {

  # Check class
  stopifnot(any(class(x) == "pwsd"))

  # Extract ACF
  acf <- x$Acf
  names <- rownames(x$BlockLength)

  # Set c to input if not NULL
  if (is.null(c)) {
    c <- x$parameters[, "c"]
  }

  # Set significance bands
  rho_k_critical <- c * sqrt(log10(x$parameters[, "n"]) / x$parameters[, "n"])

  # Plot acf
  j <- 1
  for (i in acf) {

    # Set plot title
    if (is.null(main)) {
      main = paste0("Correlogram for: ", names[j])
    }

    # Set y-axis limits
    if (is.null(ylim)) {
      ylim <- range(i$acf, rho_k_critical, -rho_k_critical)
    }

    # Plot ACF
    plot(i,
      ci = 0,
      xlab = "Lag (k)",
      main = main,
      ylim = ylim,
      ...)

    # Add implied significance bands
    graphics::abline(h = c(rho_k_critical, -rho_k_critical),
      col = "darkmagenta",
      lty = "dotdash")

    # Name Counter
    j <- j + 1

  }
}
