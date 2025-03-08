# Data setup
set.seed(32)
sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
                        n = 100, innov = rnorm(100))


# Plot HHJ ----------------------------------------------------------------
test_that("plot.hhj works", {

  x <- hhj(sim)
  expect_invisible(plot(x))
  expect_error(plot(x, iter = c("throw", 1, -1)))

  class(x) <- "testing"
  expect_error(plot(x))

})


# Plot PWSD ---------------------------------------------------------------
test_that("plot.pwsd works", {
  x <- pwsd(sim)
  expect_invisible(plot(x))
})


# Plot NPPI ---------------------------------------------------------------
test_that("plot.nppi runs without errors", {

  result <- nppi(data = sim, stat_function = mean, num_bootstrap = 100)
  expect_silent(plot.nppi(result))
})

test_that("plot.nppi handles missing JAB point values", {
  result <- list(jab_point_values = NULL)
  class(result) <- "nppi"

  expect_error(plot.nppi(result), "Error: No valid JAB point values available for plotting.")
})

test_that("plot.nppi handles all NA JAB point values", {
  result <- list(jab_point_values = rep(NA, 10))
  class(result) <- "nppi"

  expect_error(plot.nppi(result), "Error: No valid JAB point values available for plotting.")
})

test_that("plot.nppi handles missing JAB pseudo values", {
  result <- list(jab_point_values = rep(1, 10), jab_pseudo_values = NULL)
  class(result) <- "nppi"

  expect_error(plot.nppi(result), "Error: No valid JAB pseudo values available for plotting.")
})

test_that("plot.nppi handles all NA JAB pseudo values", {
  result <- list(jab_point_values = rep(1, 10), jab_pseudo_values = rep(NA, 10))
  class(result) <- "nppi"

  expect_error(plot.nppi(result), "Error: No valid JAB pseudo values available for plotting.")
})
