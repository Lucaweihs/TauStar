\name{tStar}
\alias{tStar}
\title{Computing t*}
\description{
  This functions computes the t* U-statistic for input data pairs
  (x_1,y_1), (x_2,y_2), ..., (x_n,y_n)
  using the algorithm developed by Weihs, Drton, and Leung (2015).
}
\usage{tStar(x, y, vStatistic = FALSE, resample = FALSE, numResamples = 500,
             sampleSize = min(length(x), 1000), slow = FALSE)}
\arguments{
  \item{x}{A numeric vector of x values (length >= 4).}
  \item{y}{A numeric vector of y values, should be of the same length as x.}
  \item{vStatistic}{If TRUE then will compute the V-statistic version of t*, otherwise will compute the U-Statistic version of t*. Default is to compute the U-statistic.}
  \item{resample}{If TRUE then will compute an approximation of t* using a subsettting approach: samples of size sampleSize are taken from the data numResample times, t* is computed on each subsample, and all subsample t* values are then averaged. Note that this only works when vStatistic == FALSE, in general you probably don't want to compute the V-statistic via resampling as the size of the bias depends on the sampleSize irrespective numResamples. Default is resample == FALSE so that t* is computed on all of the data, this may be slow for very large sample sizes.}
  \item{numResamples}{See resample variable description for details, this value is ignored if resample == FALSE (ignored by default).}
  \item{sampleSize}{See resample variable description for details, this value is ignored if resample == FALSE (ignored by default).}
  \item{slow}{If TRUE then will compute the t* statistic using a naive O(n^4) algorithm. This option exists to allow for comparisons of the efficient new algorithm against the old slow one. At the moment this is not implemented to work when resample == TRUE.}
}
\value{
  The numeric value of the t* statistic.
}
\references{
  Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014), no. 2, 1006--1028.
  \cr\cr
  Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
}
\examples{
\dontrun{
library(TauStar)

# Compute t* for a concordant quadruple
tStar(c(1,2,3,4), c(1,2,3,4)) # == 2/3

# Compute t* for a discordant quadruple
tStar(c(1,2,3,4), c(1,-1,1,-1)) # == -1/3

# Compute t* on random normal iid normal data
set.seed(23421)
tStar(rnorm(4000), rnorm(4000)) # near 0

# Compute t* as a v-statistic
set.seed(923)
tStar(rnorm(100), rnorm(100), vStatistic=TRUE)

# Compute an approximation of tau* via resampling
set.seed(9492)
tStar(rnorm(10000), rnorm(10000),
      resample=TRUE, sampleSize=30, numResamples=5000)
}
}
