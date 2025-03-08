# Common Data -------------------------------------------------------------
set.seed(42)
sim_data <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = 200)
small_data <- rnorm(10)  # Very small dataset


# Tests -------------------------------------------------------------------
test_that("nppi function returns expected structure", {
  result <- nppi(data = sim_data, stat_function = mean, num_bootstrap = 100)

  expect_type(result, "list")  # Ensure it's a list
  expect_s3_class(result, "nppi")  # Ensure correct S3 class

  # Check that key elements exist
  expect_named(result, c("optimal_block_length", "bias", "variance", "jab_point_values", "jab_pseudo_values", "l", "m"))

  # Check that the values are numeric
  expect_true(is.numeric(result$optimal_block_length))
  expect_true(is.numeric(result$bias))
  expect_true(is.numeric(result$variance))

  # Ensure the estimated block length is within a reasonable range
  expect_gt(result$optimal_block_length, 0)
  expect_lt(result$optimal_block_length, length(sim_data))
})

test_that("nppi stops execution when NA values are present", {
  # Copy the data and introduce an NA value
  sim_data_copy <- sim_data + 0L
  sim_data_copy[50] <- NA  # Introduce an NA value

  expect_error(
    nppi(data = sim_data_copy, stat_function = mean, num_bootstrap = 100),
    "Error: NA values detected in the time series. Please impute or remove them manually."
  )
})

test_that("nppi handles very small datasets", {
  result <- nppi(data = small_data, stat_function = mean, num_bootstrap = 100)

  expect_gt(result$optimal_block_length, 0)
  expect_lt(result$optimal_block_length, length(small_data))
})

test_that("nppi handles data.frame", {
  df <- data.frame(x = small_data)
  result <- nppi(data = df, stat_function = mean, num_bootstrap = 100)

  expect_gt(result$optimal_block_length, 0)
  expect_lt(result$optimal_block_length, length(small_data))
})

test_that("nppi handles ts object", {
  # Coerce to ts
  ts_data <- ts(small_data)
  result <- nppi(data = ts_data, stat_function = mean, num_bootstrap = 100)

  expect_gt(result$optimal_block_length, 0)
  expect_lt(result$optimal_block_length, length(small_data))
})


test_that("nppi correctly sets default values for l and m", {
  result <- nppi(data = sim_data, stat_function = mean, num_bootstrap = 100)

  expect_gt(result$l, 0)  # `l` should be a positive number
  expect_gt(result$m, 0)  # `m` should be a positive number
})

test_that("nppi jab_point_values and jab_pseudo_values are of the same length", {
  result <- nppi(data = sim_data, stat_function = mean, num_bootstrap = 100)

  expect_equal(length(result$jab_point_values), length(result$jab_pseudo_values))
})
