\name{tStar}
\alias{tStar}
\title{Computing t*}
\description{
  This functions computes the t* statistic for input data pairs
  (x_1,y_1), (x_2,y_2), ..., (x_n,y_n).
}
\usage{tStar(x, y, vStatistic = FALSE)}
\arguments{
  \item{x}{A numeric vector of x values.}
  \item{y}{A numeric vector of y values.}
  \item{vStatistic}{If TRUE then will compute the V-statistic version of t*, otherwise will compute the U-Statistic version of t*. Default is to compute the U-statistic.}
}
\value{
  The numeric value of the t* statistic.
}
\references{
  Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014), no. 2, 1006--1028.
  \cr\cr
}
\examples{
set.seed(1234)
x = rnorm(300)
y = rnorm(300)
tStar(x,y)
}
