test_that("plot.hhj works", {

  sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
    n = 50, innov = rnorm(50))

  x <- hhj(sim)
  expect_invisible(plot(x))

  expect_error(plot(x, iter = c("throw", 1, -1)))

  class(x) <- "testing"
  expect_error(plot(x))

})

test_that("plot.pwsd works", {

  sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
    n = 100, innov = rnorm(100))

  x <- pwsd(sim)

  expect_invisible(plot(x))

})
