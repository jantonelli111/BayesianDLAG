#' Plot the posterior mean and 95 percent credible interval for a main effect distributed lag
#'
#'
#' @param DLAGmodel              A fitted model from the DLAGbayes function
#' @param whichExposure          an integer indicating which exposure is to be plotted.
#'                               Must be between 1 and the number of exposures
#' @param main                   The title of the plot
#' @param ylim                   The limits of the y-axis of the plot
#' @importFrom graphics lines
#' @importFrom stats quantile
#' @export

PlotMainDLAG = function(DLAGmodel, whichExposure, main=NULL, ylim=NULL) {

  nChains = dim(DLAGmodel$betaPost)[1]
  totalScans = dim(DLAGmodel$betaPost)[2]
  dfMain = dim(DLAGmodel$betaPost)[4]
  t = DLAGmodel$t
  betaPost = DLAGmodel$betaPost
  BasisFunctions = DLAGmodel$BasisFunctions

  curve1 = array(0, dim=c(nChains, totalScans, t))

  j1 = whichExposure

  for (nc in 1 : nChains) {
    for (ni in 1 : totalScans) {
      for (k in 1 : dfMain) {
        curve1[nc,ni,] = curve1[nc,ni,] +
          betaPost[nc,ni,j1,k]*BasisFunctions[[j1]][,k]
      }

    }
  }


  if (is.null(main)) {
    main = paste("exposure", j1)
  }

  if (is.null(ylim)) {
    ylim = range(curve1)
  }

  plot(1:t, apply(curve1, 3, mean),
       ylab=expression(beta[t]), type='l', lwd=3, xlab="time",
       ylim=ylim, main=main)
  lines(1:t, apply(curve1, 3, quantile, .025), lwd=3, lty=2)
  lines(1:t, apply(curve1, 3, quantile, .975), lwd=3, lty=2)

}




#' Extract the full posterior distribution for a main effect distributed lag
#'
#'
#' @param DLAGmodel              A fitted model from the DLAGbayes function
#' @param whichExposure          an integer indicating which exposure the posterior
#'                               distribution of the main effect lag should be calculated.
#'                               Must be between 1 and the number of exposures
#'
#' @return An array containing the posterior distribution of the distributed lag. The dimensions of
#'         the array is NC x NI x T, where NC is the number of chains run, NI is the number of iterations
#'         saved per chain, and T is the number of time points we observe the exposure at.
#'
#' @export

PosteriorMainDLAG = function(DLAGmodel, whichExposure) {

  nChains = dim(DLAGmodel$betaPost)[1]
  totalScans = dim(DLAGmodel$betaPost)[2]
  dfMain = sum(!is.na(DLAGmodel$betaPost[nChains,totalScans,whichExposure,]))
  t = DLAGmodel$t
  betaPost = DLAGmodel$betaPost
  BasisFunctions = DLAGmodel$BasisFunctions

  curve1 = array(0, dim=c(nChains, totalScans, t))

  j1 = whichExposure

  for (nc in 1 : nChains) {
    for (ni in 1 : totalScans) {
      for (k in 1 : dfMain) {
        curve1[nc,ni,] = curve1[nc,ni,] +
          betaPost[nc,ni,j1,k]*BasisFunctions[[j1]][,k]
      }
    }
  }

  return(curve1)

}



#' Extract the full posterior distribution for an interaction distributed lag
#'
#'
#' @param DLAGmodel              A fitted model from the DLAGbayes function
#' @param whichExposures         a vector indicating which exposures the posterior
#'                               distribution of the interaction should be calculated.
#'                               Each must be between 1 and the number of exposures.
#'
#' @return An array containing the posterior distribution of the distributed lag. The dimensions of
#'         the array is NC x NI x T x T, where NC is the number of chains run, NI is the number of iterations
#'         saved per chain, and T is the number of time points we observe each exposure at.
#'
#' @export

PosteriorInteractionDLAG = function(DLAGmodel, whichExposures) {

  nChains = dim(DLAGmodel$betaPost)[1]
  totalScans = dim(DLAGmodel$betaPost)[2]
  dfMain = dim(DLAGmodel$betaPost)[4]
  t = DLAGmodel$t
  betaPost = DLAGmodel$betaPost
  betaIntPost = DLAGmodel$betaIntPost
  BasisFunctions = DLAGmodel$BasisFunctions
  BasisFunctionsInt = DLAGmodel$BasisFunctionsInt
  df = DLAGmodel$df
  dfI = DLAGmodel$dfI
  p = DLAGmodel$p

  intList = expand.grid(1:p, 1:p)
  w = which(intList[,1] < intList[,2])
  intList = intList[w,]

  surface = array(0, dim=c(nChains, totalScans, t, t))

  j1 = whichExposures[1]
  j2 = whichExposures[2]

  j = which(intList[,1] == j1 & intList[,2] == j2)

  for (nc in 1 : nChains) {
    for (ni in 1 : totalScans) {
      for (k in 1 : dfI[[j]]) {
        for (t1 in 1 : t) {
          for (t2 in 1 : t) {
            surface[nc,ni,t1,t2] = surface[nc,ni,t1,t2] +
              betaIntPost[nc,ni,j,k]*
              BasisFunctionsInt[[j]][t1,t2,k]
          }
        }
      }
    }
  }

  return(surface)

}
