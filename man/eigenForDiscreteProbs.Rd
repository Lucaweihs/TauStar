% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{eigenForDiscreteProbs}
\alias{eigenForDiscreteProbs}
\title{Eigenvalues for discrete asymptotic distribution}
\usage{
eigenForDiscreteProbs(p)
}
\arguments{
\item{p}{a vector of probabilities that sum to 1.}
}
\value{
the eigenvalues associated to the matrix generated by p
}
\description{
Computes the eigenvalues needed to determine the asymptotic distributions
in the mixed/discrete cases. See Nandy, Weihs, and Drton (2016)
<http://arxiv.org/abs/1602.04387> for more details.
}
