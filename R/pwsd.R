#' Politis and White (2004) Spectral Density "PWSD" Automatic Block-Length Selection
#'
#' Run the Automatic Block-Length selection method proposed by Politis and White
#'  (2004) and corrected in Patton, Politis, and White (2009). The method is
#'  based on spectral density estimation via flat-top lag windows of Politis and
#'  Romano (1995). This code was adapted from \code{\link[np]{b.star}} to add
#'  functionality and include correlogram support including an S3 method,
#'  \emph{see} Hayfield and Racine (2008).
#'
#' @param data an \eqn{n x k} data.frame, matrix, or vector (if \eqn{k = 1})
#'  where the optimal block-length will be computed for each of the \eqn{k} columns.
#' @param K_N an integer value, the maximum lags for the auto-correlation,
#'  \eqn{rho_k}, which to apply the \emph{implied hypothesis} test. Defaults to
#'  \code{max(5, log(N))}. \emph{See} Politis and White (2004) footnote c.
#' @param M_max an integer value, the upper-bound for the optimal number of lags,
#'  \eqn{M}, to compute the auto-covariance for. \emph{See} Theorem 3.3 (ii) of
#'  Politis and White (2004).
#' @param b_max a numeric value, the upper-bound for the optimal block-length.
#'  Defaults to \code{ceiling(min(3 * sqrt(n), n / 3))} per Politis and White (2004).
#' @param c a numeric value, the constant which acts as the significance level
#'  for the implied hypothesis test. Defaults to \code{qnorm(0.975)} for a 95\%
#'  confidence level. Politis and  White (2004) suggest \code{c = 2}.
#' @param round a logical value, if set to \code{FALSE} then the final block-length
#'  output will not be rounded, the default. If set to \code{TRUE} the final
#'  estimates for the optimal block-length will be rounded to whole numbers.
#' @param correlogram a logical value, if set to \code{TRUE} a plot of the
#'  correlogram (\emph{i.e.} a plot of \eqn{R(k)} vs. \eqn{k}) will be output to
#'  the console. If set to \code{FALSE}, no interim plots will be output to the
#'  console, but may be plotting later using the corresponding S3 method,
#'  \link{plot.pwsd}.
#'
#' @return an object of class 'pwsd'
#'
#' @section References:
#'
#' Andrew Patton, Dimitris N. Politis & Halbert White (2009) Correction to
#'  “Automatic Block-Length Selection for the Dependent Bootstrap” by D. Politis
#'  and H. White, Econometric Reviews, 28:4, 372-375, DOI:
#'  \doi{10.1080/07474930802459016}
#'
#' Dimitris N. Politis & Halbert White (2004) Automatic Block-Length Selection
#'  for the Dependent Bootstrap, Econometric Reviews, 23:1, 53-70, DOI:
#'  \doi{10.1081/ETC-120028836}
#'
#' Politis, D.N. and Romano, J.P. (1995), BIAS‐CORRECTED NONPARAMETRIC SPECTRAL
#'  ESTIMATION. Journal of Time Series Analysis, 16: 67-103, DOI:
#'  \doi{10.1111/j.1467-9892.1995.tb00223.x}
#'
#' Tristen Hayfield and Jeffrey S. Racine (2008). Nonparametric Econometrics:
#'  The np Package. Journal of Statistical Software 27(5). DOI:
#'  \doi{10.18637/jss.v027.i05}
#'
#' @export
#'
#' @examples
#' # Generate AR(1) time series
#' sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
#'                         n = 500, innov = rnorm(500))
#'
#' # Calculate optimal block length for series
#' pwsd(sim, round = TRUE)
#'
#' \dontrun{
#' # Use S3 Method
#' b <- pwsd(sim, round = TRUE, correlogram = FALSE)
#' plot(b)
#' }
#'
pwsd <- function(
  data,
  K_N = NULL,
  M_max= NULL,
  b_max = NULL,
  c = NULL,
  round = FALSE,
  correlogram = TRUE
  ){

  # Set PW parameters and coerce to data.frame
  data <- data.frame(data)
  n <- nrow(data)
  k <- ncol(data)

  # Save function call
  call <- match.call()

  # Set defaults pilot parameters per footnote c, page 59, Politis and White (2004)
  if (is.null(K_N)) K_N <- max(5, ceiling(log10(n)))
  if (is.null(M_max)) M_max <- ceiling(sqrt(n)) + K_N
  if (is.null(b_max)) b_max <- ceiling(min(3 * sqrt(n), n / 3)) #check PW for source
  if (is.null(c)) c <- stats::qnorm(0.975) # 95% CI

  # Store results in vectors of length k.
  b_SB_star <- rep(NA, k)
  b_CB_star <- rep(NA, k)

  # Initialize list to hold correlogram series
  autocorrelation <- vector(mode = "list", length = k)

  for(i in seq_len(length.out = k)) {

    # Construct auto-correlations, rho(k)
    autocorrelation[[i]] <- stats::acf(
      data[, i],
      lag.max = M_max,
      type = "correlation",
      plot = FALSE)

    rho.k <- autocorrelation[[i]]$acf[-1]

    # Plot correlogram
    if (isTRUE(correlogram)) {
      plot(
        autocorrelation[[i]],
        ci = stats::pnorm(c),
        ci.col =  "darkmagenta",
        xlab = "Lag (k)",
        main = paste0("Correlogram for: ", colnames(data)[i])
        )
    }

    # Set critical value for implied hypothesis test
    rho.k.critical <- c * sqrt(log10(n) / n)

    # Count number of insignificant runs for each possible k of rho(k)
    acfSignificant <- sapply(X = seq_len(length.out = (M_max - K_N + 1)),
      FUN = sigTest,
      # Parameters for sigTest
      rho.k = rho.k,
      rho.k.critical = rho.k.critical,
      K_N = K_N
      )

    # Look for insignificant runs of length K_N and take the smallest k
    if(any(acfSignificant == K_N)) {
      m_hat <- min(which(acfSignificant == K_N))
    } else {

      # If no runs of length K_N are insignificant, take smallest significant value
      if(any(abs(rho.k) > rho.k.critical)) {

        sigLags <- which(abs(rho.k) > 10)
        numSig <- length(sigLags)

        if(length(numSig) > 0) {

          m_hat <- max(sigLags)

        } else {

          # Use 1 as default if no prior conditions are satisfied
          m_hat <- 1

        }
      }
    }

    # Compute M (m_hat is at least one).
    M <- ifelse(2 * m_hat > M_max, M_max, 2 * m_hat)

    # Generate sequence of lags to look over
    k_seq <- seq(from = -M, to = M, by = 1)

    # Calculate R(k) in PW, the covariance of the series
    R.k <- stats::ccf(data[, i], data[, i],
      lag.max = M,
      type = "covariance",
      plot = FALSE)$acf

    G_hat <- sum(lambda(k_seq / M) * abs(k_seq) * R.k)

    D_CB_hat <- (4 / 3) * sum(lambda(k_seq / M) * R.k)^2
    D_SB_hat <- 2 * sum(lambda(k_seq / M) * R.k)^2

    b_SB_star[i] <- ((2 * G_hat^2) / D_SB_hat)^(1 / 3) * n^(1 / 3)
    b_CB_star[i] <- ((2 * G_hat^2) / D_CB_hat)^(1 / 3) * n^(1 / 3)

  }

  # Round final block-lengths
  if(round == FALSE) {

    b_SB_star <- ifelse(b_SB_star > b_max, b_max, b_SB_star)
    b_CB_star <- ifelse(b_CB_star > b_max, b_max, b_CB_star)

  } else {

    b_SB_star <- ifelse(b_SB_star > b_max, b_max,
      ifelse(b_SB_star < 1, 1, round(b_SB_star)))

    b_CB_star <- ifelse(b_CB_star > b_max, b_max,
      ifelse(b_CB_star < 1, 1, max(1, round(b_CB_star))))

  }

  # Collect Block-Lengths and Name rows
  blocklengths <- cbind("b_Stationary" = b_SB_star,"b_Circular" = b_CB_star)
  rownames(blocklengths) <- colnames(data)

  # Name autocorrelation list
  names(autocorrelation) <- colnames(data)

  # Collect results
  result <- structure(list(
    "BlockLength" = blocklengths,
    "Acf" = autocorrelation,
    "parameters" = cbind(n, k, c, K_N, M_max, b_max, m_hat, M, rho.k.critical),
    "Call" = call
    ), class = "pwsd")

  return(result)

}

# Helper Functions --------------------------------------------------------

# Function to construct "flat-top" lag window for spectral density estimation
## Based on Politis, D.N. and J.P. Romano (1995)
lambda <- function(t) {
  return((abs(t) >= 0) * (abs(t) < 0.5) + 2 *
      (1 - abs(t)) * (abs(t) >= 0.5) * (abs(t)<=1))
}

# Function to count the number of significant auto-correlations or given K_N
sigTest <- function(j, rho.k, rho.k.critical, K_N) {
  sum((abs(rho.k) < rho.k.critical)[j:(j + K_N - 1)])
}
