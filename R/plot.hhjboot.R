#' Plot MSE Function for HHJ Algorithm
#'
#' @param x an object of class "hhjboot"
#' @param iter a vector of \code{hhjboot()} iterations to plot. \pkg{NULL}, by
#'  default will plot all iterations
#' @param ... Arguments to be passed to methods
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' # Generate 'hhjboot' class object of optimal block length for series
#' hhj <- hhjboot(sim, sub.block.size = 10)
#'
#' ## S3 method for class 'hhjboot'
#' plot(hhj)
#' }
#'
plot.hhjboot <- function(x, iter = NULL, ...) {

  if (!inherits(x, "hhjboot")) {
    stop("use only with \"hhjboot\" objects")
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
