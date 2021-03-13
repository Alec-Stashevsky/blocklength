#' Plot MSE Function for HHJ Algorithm
#'
#' S3 Method for objects of class 'hhj'
#'
#' @param x an object of class 'hhj'
#' @param iter a vector of \code{hhj()} iterations to plot. \code{NULL}. All
#'  iterations are plotted by default.
#'
#' @inheritDotParams base::plot
#'
#'
#' @return No return value, called for side effects
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' # Generate 'hhj' class object of optimal block length for series
#' hhj <- hhj(sim, sub_sample = 10)
#'
#' ## S3 method for class 'hhj'
#' plot(hhj)
#' }
#'
plot.hhj <- function(x, iter = NULL, ...) {

  if (!inherits(x, "hhj")) {
    stop("use only with \"hhj\" objects")
  }

  p.data <- x$`MSE Data`

  if (is.null(iter)) {
    iter <- unique(p.data$Iteration)
  } else {
    stopifnot(is.numeric(iter) | is.list(iter))
  }

  # Plot all iterations
  for (i in iter) {

    # Filter on iteration
    p.data <- x$`MSE Data`[x$`MSE Data`$Iteration == i, ]

    graphics::plot.default(
      x = as.vector(p.data$BlockLength, mode = "numeric"),
      y = p.data$MSE,
      main = paste0(
        "MSE Plot for: ",
        x$Series,
        "\n",
        "Iteration: ",
        p.data$Iteration[1]
      ),
      xlab = "Block Length (l)",
      ylab = "MSE",
      col = ifelse(p.data$MSE == min(p.data$MSE), "red", "black")
    )
  }
}
