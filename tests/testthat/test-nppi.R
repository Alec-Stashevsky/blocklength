test_that("nppi function returns expected structure", {
  set.seed(42)
  sim_data <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = 200)

  result <- nppi(data = sim_data, stat_function = mean, num_bootstrap = 100)

  expect_type(result, "list")  # Ensure it's a list
  expect_s3_class(result, "nppi")  # Ensure correct S3 class

  # Check that key elements exist
  expect_named(result, c("optimal_block_length", "bias", "variance", "jab_point_values", "l", "m"))

  # Check that the values are numeric
  expect_true(is.numeric(result$optimal_block_length))
  expect_true(is.numeric(result$bias))
  expect_true(is.numeric(result$variance))

  # Ensure the estimated block length is within a reasonable range
  expect_gt(result$optimal_block_length, 0)
  expect_lt(result$optimal_block_length, length(sim_data))
})

test_that("nppi stops execution when NA values are present", {
  set.seed(42)
  sim_data <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = 200)
  sim_data[50] <- NA  # Introduce an NA value

  expect_error(
    nppi(data = sim_data, stat_function = mean, num_bootstrap = 100),
    "Error: NA values detected in the time series. Please impute or remove them manually."
  )
})

test_that("nppi handles very small datasets", {
  set.seed(42)
  small_data <- rnorm(10)  # Very small dataset

  result <- nppi(data = small_data, stat_function = mean, num_bootstrap = 100)

  expect_gt(result$optimal_block_length, 0)
  expect_lt(result$optimal_block_length, length(small_data))
})

test_that("nppi correctly sets default values for l and m", {
  set.seed(42)
  sim_data <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = 200)

  result <- nppi(data = sim_data, stat_function = mean, num_bootstrap = 100)

  expect_gt(result$l, 0)  # `l` should be a positive number
  expect_gt(result$m, 0)  # `m` should be a positive number
})
