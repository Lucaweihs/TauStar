tStar <- function(x, y, vStatistic=FALSE) {
  if(!is.numeric(x) || !is.numeric(y)) {
    stop("Input x and y to tStar must be numeric.")
  }
  if(length(x) != length(y) || length(x) < 4) {
    stop("Input x and y to tStar are of the wrong length, they must both have equal length < 4.")
  }
  if(!is.logical(vStatistic) || length(vStatistic) != 1) {
    stop("Input parameter vStatistic into function tStar must be a logical T/F value.")
  }
  ord = sort.list(x, method="quick", na.last=NA)
  x = x[ord]
  y = y[ord]
  if(vStatistic) {
    return(VTStarFastTiesRCPP(x, y))
  }
  return(TStarFastTiesRCPP(x, y))
}