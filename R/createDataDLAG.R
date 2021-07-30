#' Simulate a data set with exposures measured over time such that there is an effect
#' of exposure 1 and 2 as well as an interaction between them
#'
#'
#' @param n                   Sample size of the model
#' @param p                   Number of exposures
#' @param q                   Number of additional covariates
#' @param corr_exposures      Correlation between exposures at a point in time
#' @param corr_time           Correlation within exposures across time
#' @param beta                Optional list of true distributed lag curves for each exposure.
#'                            If left NULL, exposures 1 and 2 will have an effect on the outcome.
#' @param betaInt             Optional list of true distributed lag surfaces for each interaction
#'                            If left NULL, exposures 1 and 2 will have an interaction effect.
#' @importFrom mvtnorm stats
#' @export

createDataDLAG = function(n=200, p=5, t=37, q=5, 
                          corr_exposures = 0.5, corr_time = 0.95,
                          betaMain = NULL,
                          betaInt = NULL) {
  
  betaC = rnorm(q, sd=0.5)
  
  intList = expand.grid(1:p, 1:p)
  w = which(intList[,1] < intList[,2])
  intList = intList[w,]
  
  if (is.null(betaMain)) {
    betaMain = list()
    for (j in 1 : p) {
      betaMain[[j]] = rep(0, t)
    }
    betaMain[[1]] = (((1:t)/t)*0.25 - 0.4*((1:t)/t)^2) / 2
    
    betaMain[[2]] = c(rep(0,19), 0.00125, 0.0025, 0.005, .0075, 0.01, 0.015, 0.02, 0.03, 
                  0.04, 0.045, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0610, 0.0615, 0.0620) / 2
  } 
  
  if (is.null(betaInt)) {
    betaInt = list()
    for (j in 1 : choose(p, 2)) {
      betaInt[[j]] = matrix(0, t,t)
    }
    betaInt[[1]] = 0.7*outer(((1:t)/t)*0.6 - 0.8*((1:t)/t)^2,
                             ((1:t)/t)*0.7 - 0.7*((1:t)/t)^2) / 7
  } 
  
  residError = 1
  
  c = matrix(stats::rnorm(n*q), n, q)
  
  x = list()
  
  for (j in 1 : p) {
    x[[j]] = matrix(NA, n, t)
  }
  
  Amat = diag(corr_time, p)
  
  sigma = matrix(NA, p, p)
  
  for (rows in 1 : p) {
    for (cols in 1 : p) {
      sigma[rows,cols] = corr_exposures^(abs(rows - cols))
    }
  }
  
  allExposureData = array(NA, dim=c(n, p, t))
  allExposureData[,,1] = mvtnorm::rmvnorm(n, sigma=sigma)
  
  for (i in 1 : n) {
    for (timePoints in 2 : t) {
      allExposureData[i,,timePoints] = Amat %*% allExposureData[i,,timePoints-1] + 
        as.numeric(mvtnorm::rmvnorm(1, sigma=sigma))
    } 
  }
  
  for (j in 1 : p) {
    for (tt in 1 : t) {
      x[[j]][,tt] = (allExposureData[,j,tt] - mean(allExposureData[,j,tt]))/
        sd(allExposureData[,j,tt])
    }
  }
  
  LinPred = rep(0,n)
  for (j in 1 : p) {
    LinPred = LinPred + x[[j]] %*% betaMain[[j]]
  }
  
  LinPredInt = rep(0,n)
  for (j in 1 : choose(p,2)) {
    j1 = intList[j,1]
    j2 = intList[j,2]
    for (t1 in 1 : t) {
      for (t2 in 1 : t) {
        LinPredInt = LinPredInt + x[[j1]][,t1]*x[[j2]][,t2] * betaInt[[j]][t1,t2]
      }
    }
  }
  
  y = LinPred + LinPredInt + (c %*% betaC) + stats::rnorm(n, sd=sqrt(residError))
  
  l = list(y=y, x=x, c=c,
           betaMain=betaMain, betaInt=betaInt,
           intList)
  
  return(l)
  
}
