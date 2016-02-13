#' Efficient Computation and Testing of the t* Statistic of Bergsma and Dassios
#'
#' Computes the t* statistic corresponding to the tau star population
#' coefficient introduced by Bergsma and Dassios (Bernoulli 20(2), 2014,
#' 1006-1028) and does so in  O(n^2*log(n)) time. Also allows for
#' independence testing using the asymptotic distribution of t* as described by
#' Nandy, Weihs, and Drton (2016).
#'
#' To directly compute the t* statistic see the function tStar. If otherwise
#' interested in performing tests of independence then see the function
#' tauStarTest.
#'
#' @references
#' Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based
#' on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014), no.
#' 2, 1006--1028.
#' \cr\cr
#' Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the
#' Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
#'
#' @examples
#' \dontrun{
#' library(TauStar)
#'
#' # Compute t* for a concordant quadruple
#' tStar(c(1,2,3,4), c(1,2,3,4)) # == 2/3
#'
#' # Compute t* for a discordant quadruple
#' tStar(c(1,2,3,4), c(1,-1,1,-1)) # == -1/3
#'
#' # Compute t* on random normal iid normal data
#' set.seed(23421)
#' tStar(rnorm(4000), rnorm(4000)) # near 0
#'
#' # Compute t* as a v-statistic
#' set.seed(923)
#' tStar(rnorm(100), rnorm(100), vStatistic=TRUE)
#'
#' # Compute an approximation of tau* via resampling
#' set.seed(9492)
#' tStar(rnorm(10000), rnorm(10000),
#'       resample=TRUE, sampleSize=30, numResamples=5000)
#' }
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib TauStar
"_PACKAGE"
#> [1] "_PACKAGE"
