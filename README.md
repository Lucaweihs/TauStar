# TauStar Package

## Purpose

This package allows you to efficiently compute the U/V-statistic corresponding
to the tau* coefficient described in the paper:

Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based on a sign covariance related to Kendall's tau. Bernoulli 20 (2014), no. 2, 1006–1028.

The tau* statistic has the special property that it is 0 if and only if the
bivariate distribution it is computed upon is independent (under some weak
conditions on the bivariate distribution) and is positive otherwise. Since t*, 
the U-statistic corresponding to tau*, is an unbiased estimator of tau* this 
gives a consistent test of independence. Computing t* naively results an 
algorithm that takes O(n^4) time where n is the sample size. Luckily it is 
possible to compute t* much faster (in O(n^2*log(n)) time) using the algorithm 
described in:

Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).

This fast algorithm is implemented in this package.

## Example

A simple example of computing t* on a independent bivariate normal distribution
follows:

```
> set.seed(2342)
> x = rnorm(1000)
> y = rnorm(1000)
> tStar(x, y)
[1] 0.0003637266
```

All functionality for this package is currently wrapped in the tStar function.
