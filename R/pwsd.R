#' Politis and White (2004) Spectral Density Automatic Block-Length Selection
#'
#' Run the Automatic Block-Length selection method proposed by Politis and White
#'  (2004) and corrected in Patton, Politis, and White (2009). The method is
#'  based on spectral density estimate via flat-top lag windows of Politis and
#'  Romano (1995).
#'
#'
#' @param data an \eqn{n x k} data.frame or matrix where the optimal block-length
#'  will be computed for each of the \eqn{k} columns.
#' @param K_N the maximum lags for which to apply the \emph{implied hypothesis} testing
#'  of the autocorrelations \eqn{rho_k}. See Politis and White (2004) <INSERT>
#' @param M_max the upper bound for the optimal number of lags \eqn{M} to
#'  compute the auto-covariance. \emph{See} Theorem 3.3 (ii) of Politis and White
#' @param b_max the upper bound for the optimal block-length. Defaults to
#'  \code{ceiling(min(3 * sqrt(n), n / 3))} per Politis and White (2004).
#' @param c constant which acts as the significance level for implied
#'  hypothesis test. Defaults to \code{qnorm(0.975)} for a 95% confidence level.
#'  Politis and  White (2004) suggest \code{c=2}.
#' @param round logical. Set to \code{FALSE} by default. Setting to \code{TRUE}
#'  will round the final estimates for the optimal block-length.
#' @param correlogram logical. Setting \code{TRUE} will output a plot
#'  of the correlogram (\emph{i.e.} a plot of \eqn{R(k)} vs. \eqn{k}) to the console.
#'
#' @return an object of class 'pwsd'
#' @export
#'
#' @examples
pwsd <- function(
  data,
  K_N = NULL,
  M_max= NULL,
  b_max = NULL,
  c = NULL,
  round = FALSE,
  correlogram = TRUE
  ){

  ## Convert the data object to a data frame to handle both vectors
  ## and matrices.

  data <- data.frame(data)
  n <- nrow(data)
  k <- ncol(data)

  # Save function call
  call <- match.call()

  ## Set Defaults. Note that in footnote c, page 59, for K_N Politis
  ## and White (2004) use max(5,log10(n)). Since this must be an
  ## integer we use ceiling(log10(n)).

  if (is.null(K_N)) K_N <- max(5, ceiling(log10(n)))
  if (is.null(M_max)) M_max <- ceiling(sqrt(n)) + K_N
  if (is.null(b_max)) b_max <- ceiling(min(3 * sqrt(n), n / 3)) # not sure where this is coming from
  if (is.null(c)) c <- stats::qnorm(0.975) # 95% CI

  ## Create two vectors of length k in which we store results.

  b_SB_star <- rep(NA, k)
  b_CB_star <- rep(NA, k)

  ## Now we loop through each variable in data (i.e., column,
  ## data[,i]).

  # Initialize list to hold correlogram data
  autocorrelation <- vector(mode = "list", length = k)

  for(i in seq_len(length.out = k)) {

    ## We first obtain the autocorrelations rho(1),...,rho(M_max) (we
    ## need to drop the first autocorrelation as it is rho(0), hence
    ## acf[-1]). This is the default in acf [type="correlation"]. Note
    ## that Patton uses sample correlations after dropping the first
    ## M_max observations, while we instead use the acf to obtain
    ## rho(k).

    autocorrelation[[i]] <- stats::acf(data[, i],
      lag.max = M_max,
      type = "correlation",
      plot = FALSE)

    rho.k <- autocorrelation[[i]]$acf[-1]

    # Plot correlogram
    if (isTRUE(correlogram)) {
      plot(
        autocorrelation[[i]],
        ci = stats::pnorm(c),
        xlab = "Lag (k)",
        main = paste0("Correlogram for: ", colnames(data)[i])
        )
    }


    ## Next we compute m_hat. The use of c*sqrt(log10(n)/n) for
    ## critical values is given in footnote c of Politis and White
    ## (2004, page 59), and the approach for determining m_hat is
    ## described in footnote c.

    rho.k.critical <- c * sqrt(log10(n) / n)

    ## Compute the number of insignificant runs following each rho(k),
    ## k=1,...,M_max.

    acfSignificant <- sapply(X = seq_len(length.out = (M_max - K_N + 1)),
      FUN = sigTest,
      # Parameters for sigTest
      rho.k = rho.k,
      rho.k.critical = rho.k.critical,
      K_N = K_N
      )

    ## If there are any values of rho(k) for which the K_N proceeding
    ## values of rho(k+j), j=1,...,K_N are all insignificant, take the
    ## smallest rho(k) such that this holds (see footnote c for
    ## further details).

    if(any(acfSignificant == K_N)) {
      m_hat <- min(which(acfSignificant == K_N))
    } else {

      ## If no runs of length K_N are insignificant, take the smallest
      ## value of rho(k) that is significant.

      if(any(abs(rho.k) > rho.k.critical)) {

        sigLags <- which(abs(rho.k) > 10)
        numSig <- length(sigLags)

        if(length(numSig) > 0) {
          m_hat <- max(sigLags)
        } else {
          m_hat <- 1
        }

# OLD ---------------------------------------------------------------------


        #if(numSig == 1) {

          ## When only one lag is significant, m_hat is the sole
          ## significant rho(k).

          #m_hat <- sigLags

      #  } else {

          ## If there are more than one significant lags but no runs
          ## of length K_N, take the largest value of rho(k) that is
          ## significant.

       #   m_hat <- max(sigLags)

      #  }

    #  } else {

        ## When there are no significant lags, m_hat must be the
        ## smallest positive integer (footnote c), hence m_hat is set
        ## to one.

     #   m_hat <- 1



      }

    }

    ## Compute M (m_hat is at least one).

    M <- ifelse(2 * m_hat > M_max, M_max, 2 * m_hat)

    ## We compute b_SB_star and b_CB_star using the formulas in the above
    ## references. Now we require the autocovariance R(k) (hence
    ## type="covariance" in the acf call). Note that Patton uses
    ## sample covariances after dropping the first M_max observations,
    ## while we instead use the acf with type="covariance" to obtain
    ## R(k). Note also that we require R(0) hence we do not drop it as
    ## we did for rho(k) via acf(...)$acf[-1].

    k_seq <- seq(from = -M, to = M, by = 1)

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

  ## The user can choose whether they want rounded values returned or
  ## not. b_CB_star is rounded up, b_SB_star simply rounded but both must
  ## be positive integers.

  if(round == FALSE) {

    b_SB_star <- ifelse(b_SB_star > b_max, b_max, b_SB_star)
    b_CB_star <- ifelse(b_CB_star > b_max, b_max, b_CB_star)

  } else {

    b_SB_star <- ifelse(b_SB_star > b_max, b_max,
      ifelse(b_SB_star < 1, 1, round(b_SB_star)))

    b_CB_star <- ifelse(b_CB_star > b_max, b_max,
      ifelse(b_CB_star < 1, 1, max(1,round(b_CB_star))))

  }

  # Collect results
  result <- structure(list(
    "Block-Length" = cbind("b_Stationary" = b_SB_star,"b_Circular" = b_CB_star),
    "Acf" = autocorrelation,
    "Call" = call
    ), class = c("pwsd", "list"))

  return(result)

}


## $Id: ppw.R,v 1.47 2008/12/12 14:52:17 jracine Exp jracine $

## Original code in Matlab by A. Patton, R translation and
## modifications by C. Parmeter and J. Racine.
##
## We are grateful to Andrew Patton and Dimitris Politis for their
## assistance and feedback. Kindly report features, deficiencies, and
## improvements to racinej@mcmaster.ca.
##
## The citation is A. Patton, D.N. Politis, and H. White (2008,
## forthcoming), "CORRECTION TO `Automatic Block-Length Selection for
## the Dependent Bootstrap' by D.N. Politis and H. White". This is
## based on the article by Politis, D.N., and H. White (2004),
## "Automatic block-length selection for the dependent bootstrap."
## Econometric Reviews, vol. 23.
##
## INPUTS:  data, an n x k matrix.
##
## OUTPUTS: b.star, a 2 x k vector of optimal bootstrap block lengths
## for the stationary bootstrap and circular bootstrap (b_SB_star,
## b_CB_star).

## The function lambda() is used to construct a "flat-top" lag window for
## spectral estimation based on Politis, D.N. and J.P. Romano (1995),
## "Bias-Corrected Nonparametric Spectral Estimation", Journal of Time
## Series Analysis, vol. 16, No. 1.


# Helper Functions --------------------------------------------------------

# “flat-top” lag-window of Politis and Romano (1995) - Trapazoidal Shape symmetric about zero.
lambda <- function(t) {
  return((abs(t) >= 0) * (abs(t) < 0.5) + 2 *
      (1 - abs(t)) * (abs(t) >= 0.5) * (abs(t)<=1))
}

sigTest <- function(j, rho.k, rho.k.critical, K_N) {
  sum((abs(rho.k) < rho.k.critical)[j:(j + K_N - 1)])
}


## The function b.star() returns the optimal bootstrap block
## lengths. Note that an example for usage appears at the bottom of
## this file. If you use this function as input into a routine such as
## tsboot() in the boot library (Angelo Canty and Brian Ripley
## (2008). boot: Bootstrap R (S-Plus) Functions. R package version
## 1.2-34.) you ought to use the option round=TRUE.



# ALEC NOTES:

# M is the lag-length for the autocorrelation function ,rho.k, which is used to estimate the m_hat = M / 2. An 'implied significance' test is used to find the significant autocorrelations up to m_max (set by Politis formulas)


# rho.k is the autocorrelation and R.k is the autocovariance!

# g_hat is the spectral density function
