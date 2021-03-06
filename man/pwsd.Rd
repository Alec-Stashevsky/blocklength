% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwsd.R
\name{pwsd}
\alias{pwsd}
\title{Politis and White (2004) Spectral Density "PWSD" Automatic Block-Length Selection}
\usage{
pwsd(
  data,
  K_N = NULL,
  M_max = NULL,
  m_hat = NULL,
  b_max = NULL,
  c = NULL,
  round = FALSE,
  correlogram = TRUE
)
}
\arguments{
\item{data}{an \eqn{n x k} data.frame, matrix, or vector (if \eqn{k = 1})
where the optimal block-length will be computed for each of the \eqn{k}
columns.}

\item{K_N}{an integer value, the maximum lags for the auto-correlation,
\eqn{rho_k}, which to apply the \emph{implied hypothesis} test. Defaults to
\code{max(5, log(N))}. \emph{See} Politis and White (2004) footnote c.}

\item{M_max}{an integer value, the upper-bound for the optimal number of lags,
\eqn{M}, to compute the auto-covariance for. \emph{See} Theorem 3.3 (ii) of
Politis and White (2004).}

\item{m_hat}{an integer value, if set to \code{NULL} (the default), then
\code{m_hat} is estimated as the smallest integer after which the correlogram
appears negligible for \code{K_N} lags. In problematic cases, setting
\code{m_hat} to an integer value can be used to override the estimation
procedure.}

\item{b_max}{a numeric value, the upper-bound for the optimal block-length.
Defaults to \code{ceiling(min(3 * sqrt(n), n / 3))} per Politis and White
(2004).}

\item{c}{a numeric value, the constant which acts as the significance level
for the implied hypothesis test. Defaults to \code{qnorm(0.975)} for a
two-tailed 95\% confidence level. Politis and  White (2004) suggest
\code{c = 2}.}

\item{round}{a logical value, if set to \code{FALSE} then the final
block-length output will not be rounded, the default. If set to \code{TRUE}
the final estimates for the optimal block-length will be rounded to whole
numbers.}

\item{correlogram}{a logical value, if set to \code{TRUE} a plot of the
correlogram (\emph{i.e.} a plot of \eqn{R(k)} vs. \eqn{k}) will be output to
the console. If set to \code{FALSE}, no interim plots will be output to the
console, but may be plotted later using the corresponding S3 method,
\link{plot.pwsd}.}
}
\value{
an object of class 'pwsd'
}
\description{
Run the Automatic Block-Length selection method proposed by Politis and White
 (2004) and corrected in Patton, Politis, and White (2009). The method is
 based on spectral density estimation via flat-top lag windows of Politis and
 Romano (1995). This code was adapted from \code{\link[np]{b.star}} to add
 functionality and include correlogram support including an S3 method,
 \emph{see} Hayfield and Racine (2008).
}
\section{References}{


Andrew Patton, Dimitris N. Politis & Halbert White (2009) Correction to
     "Automatic Block-Length Selection for the Dependent Bootstrap" by D.
     Politis and H. White, Econometric Review, 28:4, 372-375, DOI:
     \doi{10.1080/07474930802459016}

Dimitris N. Politis & Halbert White (2004) Automatic Block-Length Selection
     for the Dependent Bootstrap, Econometric Reviews, 23:1, 53-70, DOI:
     \doi{10.1081/ETC-120028836}

Politis, D.N. and Romano, J.P. (1995), Bias-Corrected Nonparametric Spectral
     Estimation. Journal of Time Series Analysis, 16: 67-103, DOI:
     \doi{10.1111/j.1467-9892.1995.tb00223.x}

Tristen Hayfield and Jeffrey S. Racine (2008). Nonparametric Econometrics:
     The np Package. Journal of Statistical Software 27(5). DOI:
     \doi{10.18637/jss.v027.i05}
}

\examples{
# Generate AR(1) time series
sim <- stats::arima.sim(list(order = c(1, 0, 0), ar = 0.5),
                        n = 500, innov = rnorm(500))

# Calculate optimal block length for series
pwsd(sim, round = TRUE)


# Use S3 Method
b <- pwsd(sim, round = TRUE, correlogram = FALSE)
plot(b)

}
