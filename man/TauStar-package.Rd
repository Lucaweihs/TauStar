\name{TauStar-package}
\alias{TauStar-package}
\alias{TauStar}
\docType{package}
\title{
Efficient Computation of the t* Statistic of Bergsma and Dassios
}
\description{
Computes the t* statistic corresponding to the tau* population coefficient introduced by Bergsma and Dassios (Bernoulli 20(2), 2014, 1006-1028) and does so in O(n^2*log(n)) time using the algorithm described by Weihs, Drton, and Leung (2015). Can provide both the V-statistic and U-statistic related to the tau star measure depending on user preference.
}
\details{
\tabular{ll}{
Package: \tab TauStar\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-30\cr
License: \tab GPL (>= 2)\cr
}
% An overview of how to use the package, including the most important functions
All functionality for this package is accessed through the tStar function.
}
\author{
Luca Weihs

Maintainer: Luca Weihs <lucaw@uw.edu>
}
\references{
Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014), no. 2, 1006--1028.
\cr\cr
Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
}
\keyword{ package }
%\seealso{
% Optional links to other man pages, e.g.
% \code{\link[<pkg>:<pkg>-package]{<pkg>}}
%}
\examples{
% simple examples of the most important functions
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
