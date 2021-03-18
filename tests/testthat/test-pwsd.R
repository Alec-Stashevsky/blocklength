test_that("flat-top lag window (lambda function) works", {

  t <- runif(10, min = 0.5000001, max = 1)
  t1 <- runif(5, min = 0, max = 0.5)

  # Edge Cases
  expect_equal(lambda(0), 1)
  expect_equal(lambda(0.5), 1)
  expect_equal(lambda(2), 0)
  expect_equal(lambda(NA), NA_integer_)

  # Vector Inputs
  expect_equal(lambda(t), 2 * (1 - abs(t)))
  expect_equal(lambda(t1), rep(1, length(t1)))

})

test_that("implied hypothesis test works", {

  sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
    n = 100, innov = rnorm(100))

  n <- length(sim)

  autocorrelation <- stats::acf(
    sim,
    lag.max = 20,
    type = "correlation",
    plot = FALSE)

  # Remove lag 0 for testing
  rho_k <- autocorrelation$acf[-1]
  rho_k_critical <- 2 * sqrt(log10(n) / n)

  expect_equal(
    implied_hypothesis(rho_k = rho_k, rho_k_critical = 0, K_N = 5),
    length(rho_k)
    )

  expect_equal(
    implied_hypothesis(rho_k = rho_k, rho_k_critical = 1, K_N = 5), 1
    )

})


test_that("pwsd works", {

  sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
    n = 100, innov = rnorm(100))

  expect_type(pwsd(sim), "list")
  expect_is(pwsd(sim), "pwsd")
  expect_error(pwsd(NA))

})
