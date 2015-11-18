# Copyright (C) 2015 Luca Weihs
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

tStar <- function(x, y, vStatistic = FALSE, resample = FALSE,
                  numResamples = 500, sampleSize = min(length(x), 1000),
                  slow = FALSE) {
  if(!is.numeric(x) || !is.numeric(y)) {
    stop("Input x and y to tStar must be numeric.")
  }
  if(length(x) != length(y) || length(x) < 4) {
    stop("Input x and y to tStar are of the wrong length, they must both have equal length < 4.")
  }
  if(!is.logical(vStatistic) || length(vStatistic) != 1) {
    stop("Input parameter vStatistic into function tStar must be a logical T/F value.")
  }
  if(!is.logical(slow) || length(slow) != 1) {
    stop("Input parameter slow into function tStar must be a logical T/F value.")
  }
  if(!is.logical(resample) || length(resample) != 1) {
    stop("Input parameter resample into function tStar must be a logical T/F value.")
  }
  if(resample && numResamples <= 0) {
    stop("When resampling the number of resamples must be positive.")
  }
  if(resample && (sampleSize < 4 || sampleSize > length(x))) {
    stop("When resampling the sample size must be greater than 3 and less than the length of the input data.")
  }
  if(resample && slow) {
    stop("Resampling is not currently implemented with the slow algorithm.")
  }
  if(resample && vStatistic) {
    stop("Resampling is not currently implemented when computing V-statistics. Note that you probably don't want to compute the V-statistic via resampling as the size of the bias would depend on the size of subsets chosen independent of the number of resamples.")
  }
  ord = sort.list(x, method="quick", na.last=NA)
  x = x[ord]
  y = y[ord]
  if(resample) {
    return(TStarFastResampleRCPP(x, y, numResamples, sampleSize))
  } else if(slow) {
    return(TStarSlowTiesRCPP(x, y, vStatistic))
  } else if (vStatistic) {
    return(VTStarFastTiesRCPP(x, y))
  } else {
    return(TStarFastTiesRCPP(x, y))
  }
}

qHoeffInd <- function(x) {
  if (x == 0) {
    return(-1)
  }
  if (x == 1) {
    return(Inf)
  }
  binarySearch = function(x, lastLeft, lastRight, error = 10^-4) {
    mid = (lastLeft + lastRight) / 2
    pMid = pHoeffInd(mid)
    if (lastRight - lastLeft < error) {
      return(mid)
    } else if(pMid > x) {
      return(binarySearch(x, lastLeft, mid, error))
    } else {
     return(binarySearch(x, mid, lastRight, error))
    }
  }
  right = 1
  while(pHoeffInd(right) < x) {
    right = right * 2
  }
  return(binarySearch(x, -1, right))
}

rHoeffInd <- function(n) {
  sims = numeric(n)
  for(i in 1:50) {
    for(j in 1:50) {
      sims = sims + (36 / pi^4) * (1 / (i^2 * j^2)) * (rchisq(n, df=1) - 1)
    }
  }
  return(sims)
}

pHoeffInd <- function(x, error = 10^-5) {
  if (length(x) != 1) {
    return(sapply(x, function(y) { HoeffIndCdfRCPP(y + 1, error) }))
  }
  return(HoeffIndCdfRCPP(x + 1, error))
}

dHoeffInd <- function(x, error = 1/2 * 10^-3) {
  if (length(x) != 1) {
    return(sapply(x, function(y) { HoeffIndPdfRCPP(y + 1, error) }))
  }
  HoeffIndPdfRCPP(x + 1, error)
}

pDisHoeffInd <- function(x, p, q, error = 10^-5) {
  eigenP = eigenForDiscreteProbs(p)
  eigenQ = eigenForDiscreteProbs(q)
  return(HoeffIndDiscreteCdfRCPP(x + 4 * sum(eigenP) * sum(eigenQ),
                                 eigenP,
                                 eigenQ,
                                 error))
}

dDisHoeffInd <- function(x, p, q, error = 10^-3) {
  eigenP = eigenForDiscreteProbs(p)
  eigenQ = eigenForDiscreteProbs(q)
  return(HoeffIndDiscretePdfRCPP(x + 4 * sum(eigenP) * sum(eigenQ),
                                 eigenP,
                                 eigenQ,
                                 error))
}

rDisHoeffInd <- function(n, p, q) {
  eigenP = eigenForDiscreteProbs(p)
  eigenQ = eigenForDiscreteProbs(q)
  asymResults = numeric(asymSims)
  for(i in 1:length(eigenP)) {
    for(j in 1:length(eigenQ)) {
      asymResults = asymResults + 4 * eigenP[i] * eigenQ[j] * (rchisq(n, df=1) - 1)
    }
  }
  return(asymResults)
}

print.tstest <- function(tsObj) {
  if (tsObj$mode %in% c("continuous", "discrete")) {
    cat(paste("Test Type: asymptotic", tsObj$mode, "\n"))
  } else {
    cat(paste("Test Type: permutation (", tsObj$resamples," simulations)\n", sep=""))
  }
  cat(paste("Input Length:", length(tsObj$x), "\n\n"))

  cat(paste("Results:\n"))
  df = data.frame(round(tsObj$tStar, 5))
  if (tsObj$mode %in% c("continuous", "discrete")) {
    df = cbind(df, round(tsObj$pVal, 5))
    colnames(df) = c("t* value", "Asym. p-val")
  } else {
    df = cbind(df, round(tsObj$permPVal, 5))
    colnames(df) = c("t* value", "Perm. p-val")
  }
  row.names(df) = ""
  print(df)
}

tauStarTest <- function(x, y, mode="continuous", resamples = 1000) {
  uniqueX = unique(x)
  uniqueY = unique(y)

  toReturn = list()
  class(toReturn) = "tstest"
  toReturn$mode = mode
  toReturn$x = x
  toReturn$y = y
  toReturn$tStar = tStar(x, y)
  toReturn$resamples = resamples

  n = length(x)
  if (mode == "continuous") {
    if (length(uniqueX) != length(x) || length(uniqueY) != length(y)) {
      stop("Input vectors to tauStarTest have duplicates (repeated entries) while the mode is set to continuous.")
    }
    toReturn$pVal = 1 - pHoeffInd(n * toReturn$tStar)

  } else if (mode == "discrete") {
    p = as.numeric(table(x)) / n
    q = as.numeric(table(y)) / n
    toReturn$pVal = 1 - pDisHoeffInd(n * toReturn$tStar, p=p, q=q)

  } else if (mode == "permutation") {
    sampleTStars = numeric(resamples)
    for (i in 1:resamples) {
      sampleTStars[i] = tStar(sample(x), y)
    }
    toReturn$permPVal = mean(sampleTStars >= toReturn$tStar)
  } else {
    stop("Invalid mode as input to tauStarTest.")
  }
  return(toReturn)
}




