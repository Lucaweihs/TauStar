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
