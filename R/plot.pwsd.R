#' Plot Correlogram for Politis and White Auto\eqn{-}Correlation Implied Hypothesis Test
#'
#' S3 Method for objects of class 'pwsd'
#' \emph{See} \code{?plot.acf} of the \pkg{stats} package for more customization
#'  options on the correlogram, from which \code{plot.pwsd} is based
#'
#' @param x an of object of class 'pwsd' or 'acf'
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
plot.pwsd <- function(x, c = NULL, main = NULL, ...) {

  # Check class
  stopifnot(any(class(x) == "pwsd"))

  # Extract ACF
  acf <- x$Acf
  names <- rownames(x$BlockLength)

  # Set c to input if not NULL
  if (is.null(c)) {
    c <- x$parameters[, "c"]
  } else {c <- c}

  # Plot acf
  j <- 1
  for (i in acf) {

    if (is.null(main)) {
      main = paste0("Correlogram for: ", names[j])
    } else {
      main = main
    }

    # Plot ACT
    plot(i,
      ci = stats::pnorm(c),
      ci.col =  "darkmagenta",
      xlab = "Lag (k)",
      main = main,
      ...)

    # Name Counter
    j <- j + 1

  }
}
