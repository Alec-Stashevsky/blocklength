test_that("hhjMSE function works", {

  series <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
    n = 100, innov = rnorm(100))

  n <- length(series)

  m <- round(n^(1 / 5) * n^(1 / 3))

  # Initialize overlapping sub samples list
  series.list <- vector(mode = "list", length = length((n - m + 1)))

  for (i in seq_len(length.out = (n - m + 1))) {
    series.list[[i]] <- series[seq(from = i, to = (i + m - 1), by = 1)]
  }

  # Tests
  expect_gte(
    hhjMSE(1, n = n, m = m,
      series.list = series.list, nb = 5,
      bofb = 1L, v_star = 2),
    0)

  expect_error(
    hhjMSE(0, n = n, m = m,
      series.list = series.list, nb = 5,
      bofb = 1L, v_star = 2)
    )

  expect_error(
    hhjMSE(-1, n = n, m = m,
      series.list = series.list, nb = 5,
      bofb = 1L, v_star = 2)
  )

  expect_equal(
    hhjMSE(1, n = n, m = m,
      series.list = series, nb = 5,
      bofb = 1L, v_star = 2),
    NA_integer_)

})


test_that("hhj function works", {

  series <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
    n = 20, innov = rnorm(20))

  expect_type(hhj(series), "list")
  expect_type(
    suppressWarnings(
      hhj(series, sub_sample = 2, k = "one-sided")), "list"
    )
  expect_type(
    suppressWarnings(
      hhj(series, search_grid = 5, k = "bias/variance")), "list"
    )
  expect_is(hhj(series, grid_step = 4, pilot_block_length = 3), "hhj")
  expect_error(hhj(NA))
  expect_message(hhj(series, plots = FALSE))

  expect_error(hhj(series, grid_step = 3.3))
  expect_error(hhj(series, grid_step = c(3, 2, 1)))
  expect_error(hhj(series, search_grid = 3.3))
  expect_error(hhj(series, sub_sample = length(series) + 10))
  expect_error(hhj(series, k = "testing"))

})

test_that("hhj parrallelization works", {

  series <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
    n = 20, innov = rnorm(20))

  cl <- parallel::makeCluster(2)

  expect_is(hhj(series, cl = cl), "hhj")

  parallel::stopCluster(cl)

})
