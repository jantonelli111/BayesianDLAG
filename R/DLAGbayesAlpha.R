#' Estimate effects of repeated exposures using distributed lag models and controls average posterior
#' inclusion probabilities for null associations
#'
#'
#' @param y                   The outcome to be analyzed
#' @param x                   A list where each element is a matrix containing the values
#'                            of each exposure over time
#' @param c                   An n by q matrix of additional, time invariant covariates to adjust for
#' @param nScans              The number of MCMC scans to run
#' @param nBurn               The number of MCMC scans that will be dropped as a burn-in
#' @param thin                This number represents how many iterations between each scan
#'                            that is kept
#' @param PCAthresholdBasis   This number represents what percentage of variation we want our basis functions
#'                            from PCA to capture from the original exposures
#' @param PCAthresholdMain    The proportion of variation to keep in the main effect DLAG curves. We recommend
#'                            leaving this parameter fixed at 1
#' @param PCAthresholdInt     The proportion of variation to keep in the interaction DLAG surfaces We recommend
#'                            leaving this parameter fixed at 0.999
#' @param VariableSel         Whether to perform variable selection on main effect distributed lag surfaces
#' @param Basis               Indicates what type of basis functions to use. The default is PCA, which then
#'                            uses the PCA threshold variable to determine how many basis functions to keep.
#'                            The user may also input "Cubic" to get Cubic polynomial basis functions
#' @param alphaMain           The level at which we aim to restrict the posterior inclusion probabilities
#'                            for main effects with no association with the outcome
#' @param alphaInt            The level at which we aim to restrict the posterior inclusion probabilities
#'                            for interactions with no association with the outcome
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom stats lm rgamma runif rbeta cov var rnorm coef
#' @return A list of values that contain the posterior inclusion probabilities for the main effect
#'         distributed lags and interaction distributed lag surfaces. The list also contains the
#'         posterior distribution for the distributed lag parameters as well as the basis functions
#'         used to estimate their distributed lag function. These posterior distributions do not
#'         represent the distributed lag surfaces themselves, but these can be obtained using the
#'         other functions in the package. See PosteriorMainDLAG() and PosteriorInteractionDLAG().
#'
#' @export



DLAGbayesAlpha = function(y, x, c=NULL,
                              nScans = 2000, nBurn = 1000, thin = 2,
                              PCAthresholdBasis = 0.95,
                              PCAthresholdMain = 1,
                              PCAthresholdInt = 0.999,
                              VariableSel = TRUE, Basis = "PCA",
                              alphaMain = 0.05, alphaInt = 0.01) {

  print("Setting up MCMC")

  n = length(y)
  p = length(x)
  t = dim(x[[1]])[2]
  q = 0
  designC = cbind(rep(1,n))
  if (is.null(c) == FALSE) {
    c = as.matrix(c)
    q = dim(c)[2]
    
    designC = cbind(rep(1,n), as.matrix(c))
  }
  CC <- t(designC) %*% designC
  
  SigmaC = diag(10000, q+1)
  SigmaCInv <- chol2inv(chol(SigmaC))
  
  intList = expand.grid(1:p, 1:p)
  w = which(intList[,1] < intList[,2])
  intList = intList[w,]
  
  ## Create basis functions for main effects and interactions
  if (Basis == "PCA") {
    createMats = CreatePCA(x=x, PCAthresholdBasis = PCAthresholdBasis)
  } else if (Basis == "Cubic") {
    createMats = CreateCubic(x=x)
  } else {
    createMats = CreateCombinedCubic(x=x, PCAthresholdBasis = PCAthresholdBasis)
  }
  
  designMain = createMats$designMain
  designInt = createMats$designInt
  BasisFunctions = createMats$BasisFunctions
  BasisFunctionsInt = createMats$BasisFunctionsInt
  df = createMats$df
  dfI = createMats$dfI
  dfMain = createMats$dfMain
  dfInt = createMats$dfInt
  
  overallDesign = designMain[,1:(df[[1]] + 1),1]
  for (j in 2 : p) {
    overallDesign = cbind(overallDesign, designMain[,1:(df[[j]] + 1),j])
  }
  for (j in 1 : choose(p, 2)) {
    overallDesign = cbind(overallDesign,designInt[,1 : dfI[[j]],j])
  }
  
  overallDesign = scale(overallDesign)
  
  if (is.null(c) == TRUE) {
    if (n > (2*ncol(overallDesign))) {
      ggModel = lm(y ~ overallDesign)
      sigmaSqEst = mean((y - cbind(rep(1, n), overallDesign) %*%
                           ggModel$coefficients)^2)
      
      y = y/sqrt(sigmaSqEst)
      designMain = designMain/sqrt(sigmaSqEst)
      designInt = designInt/sqrt(sigmaSqEst)
    } else {
      ggModel = glmnet::cv.glmnet(y=y, x=overallDesign)
      sigmaSqEst = mean((y - cbind(rep(1, n), overallDesign) %*%
                           as.numeric(coef(ggModel, s="lambda.min")))^2)
      
      y = y/sqrt(sigmaSqEst)
      designMain = designMain/sqrt(sigmaSqEst)
      designInt = designInt/sqrt(sigmaSqEst)
    }
  } else {
    
    
    if (n > (2*ncol(overallDesign))) {
      ggModel = lm(y ~ overallDesign + as.matrix(c))
      sigmaSqEst = mean((y - cbind(rep(1, n), overallDesign, as.matrix(c)) %*%
                           ggModel$coefficients)^2)
      
      y = y/sqrt(sigmaSqEst)
      designMain = designMain/sqrt(sigmaSqEst)
      designInt = designInt/sqrt(sigmaSqEst)
      c = as.matrix(c)/sqrt(sigmaSqEst)
    } else {
      ggModel = glmnet::cv.glmnet(y=y, x=cbind(overallDesign, as.matrix(c)))
      sigmaSqEst = mean((y - cbind(rep(1, n), overallDesign, as.matrix(c)) %*%
                           as.numeric(coef(ggModel, s="lambda.min")))^2)
      
      y = y/sqrt(sigmaSqEst)
      designMain = designMain/sqrt(sigmaSqEst)
      designInt = designInt/sqrt(sigmaSqEst)
      c = as.matrix(c)/sqrt(sigmaSqEst)
    }
  }
  
  ## Scale the data appropriately
  meanVecMain = apply(designMain, 2:3, mean, na.rm=TRUE)
  meanVecInt = apply(designInt, 2:3, mean, na.rm=TRUE)
  
  sdVecMain = apply(designMain, 2:3, sd, na.rm=TRUE)
  sdVecInt = apply(designInt, 2:3, sd, na.rm=TRUE)
  
  for (j in 1 : p) {
    for (k in 1 : (df[[j]] + 1)) {
      designMain[,k,j] = (designMain[,k,j] - meanVecMain[k,j]) / sdVecMain[k,j]
    }
  }
  
  for (j in 1 : choose(p, 2)) {
    for (k in 1 : dfI[[j]]) {
      designInt[,k,j] = (designInt[,k,j] - meanVecInt[k,j]) / sdVecInt[k,j]
    }
  }
  
  QmatMain = DvecMain = df2 = list()
  QmatMain2 = DvecMain2 = list()
  designMain2 = array(NA, dim=dim(designMain))
  for (j in 1 : p) {
    tempX1 = designMain[,1,j]
    tempX2 = designMain[,2 : (df[[j]] + 1),j]
    SVD = svd((1/n)*t(tempX2) %*% tempX2)
    nComponents = which((cumsum(SVD$d) / sum(SVD$d)) >= PCAthresholdMain)[1]
    df2[[j]] = nComponents + 1
    QmatMain[[j]] = as.matrix(SVD$u[,1:nComponents])
    DvecMain[[j]] = SVD$d[1:nComponents]
    
    if (nComponents == 1) {
      tempDesign = cbind(tempX1, tempX2 %*% QmatMain[[j]] %*% (1/sqrt(DvecMain[[j]])))
    } else {
      tempDesign = cbind(tempX1, tempX2 %*% QmatMain[[j]] %*% diag(1/sqrt(DvecMain[[j]])))
    }
    
    SVD2 = svd((1/n)*t(tempDesign) %*% tempDesign)
    QmatMain2[[j]] = as.matrix(SVD2$u)
    DvecMain2[[j]] = SVD2$d
    
    designMain2[,1:df2[[j]],j] = tempDesign %*% QmatMain2[[j]] %*% diag(1/sqrt(DvecMain2[[j]])) 
  }
  
  QmatInt = DvecInt = dfI2 = list()
  QmatInt2 = DvecInt2 = list()
  designInt2 = array(NA, dim=dim(designInt))
  for (j in 1 : choose(p,2)) {
    tempX1 = designInt[,1,j]
    tempX2 = designInt[,2:dfI[[j]],j]
    SVD = svd((1/n)*t(tempX2) %*% tempX2)
    nComponents = which((cumsum(SVD$d) / sum(SVD$d)) >= PCAthresholdInt)[1]
    dfI2[[j]] = nComponents+1
    QmatInt[[j]] = as.matrix(SVD$u[,1:nComponents])
    DvecInt[[j]] = SVD$d[1:nComponents]
    
    designInt2[,1,j] = tempX1
    
    if (nComponents == 1) {
      tempDesign = cbind(tempX1, tempX2 %*% QmatInt[[j]] %*% (1/sqrt(DvecInt[[j]])))
    } else {
      tempDesign = cbind(tempX1, tempX2 %*% QmatInt[[j]] %*% diag(1/sqrt(DvecInt[[j]])))
    }
    
    SVD2 = svd((1/n)*t(tempDesign) %*% tempDesign)
    QmatInt2[[j]] = as.matrix(SVD2$u)
    DvecInt2[[j]] = SVD2$d
    
    designInt2[,1:dfI2[[j]],j] = tempDesign %*% QmatInt2[[j]] %*% diag(1/sqrt(DvecInt2[[j]]))
  }
  
  nChains = 2
  
  gammaPost = array(NA, dim=c(nChains, nScans, p))
  gammaIntPost = array(NA, dim=c(nChains, nScans, choose(p,2)))
  betaPost = array(NA, dim=c(nChains, nScans, p, dfMain))
  betaIntPost = array(NA, dim=c(nChains, nScans, choose(p,2), dfInt))
  betaCPost = array(NA, dim=c(nChains, nScans, q+1))
  sigmaPost = array(NA, dim=c(nChains, nScans))
  sigmaBetaPost = array(NA, dim=c(nChains, nScans))
  sigmaBetaIntPost = array(NA, dim=c(nChains, nScans))
  
  ## TODO MAKE THESE RANDOM
  gammaPost[,1,] = 0
  gammaIntPost[,1,] = 0
  betaPost[,1,,] = 0
  betaCPost[,1,] = 0
  betaIntPost[,1,,] = 0
  sigmaPost[,1] = rgamma(nChains, 1, 1)
  sigmaBetaPost[,1] = rgamma(nChains, 1, 1)
  sigmaBetaIntPost[,1] = rgamma(nChains, 1, 1)
  
  aSig = 0.001
  bSig = 0.001
  
  aTau = 3
  bTau = p
  
  aSigBeta = 0.5
  bSigBeta = 0.5
  
  TempLinPred = array(0, dim=c(nChains, n, p))
  TempLinPredInt = array(0, dim=c(nChains, n, choose(p,2)))
  
  logTauVec = c(seq(-20, -10, by=1),
                seq(-9, -1, by=0.2),
                seq(-0.9, -0.1, by=0.1),
                seq(-0.1, -0.01, by=0.01),
                seq(-0.01, -0.001, by=0.001),
                seq(-0.001, -0.0001, by=0.0001),
                seq(-0.0001, -0.00001, by=0.00001),
                seq(-0.00001, -0.000001, by=0.000001),
                seq(-0.000001, -0.0000001, by=0.0000001),
                seq(-0.0000001, -0.00000001, by=0.00000001),
                seq(-0.00000001, -0.000000001, by=0.000000001),
                seq(-0.000000001, -0.0000000001, by=0.0000000001),
                seq(-0.0000000001, -0.00000000001, by=0.00000000001))
  tauVec = exp(logTauVec)
  
  ## Finding tau as a function of sigma values for main effects
  sigmaBetaVec = c(0.0000001, seq(0.0001, 0.001, length=5),
                   seq(0.002, 0.01, length=5),
                   seq(0.02, 0.1, length=5),
                   seq(0.2, 1, length=5),
                   3,5,10,100,10000)
  sigmaVec = c(seq(as.numeric(0.01*var(y)), as.numeric(var(y)), length=20),
               as.numeric(1.5*var(y)),
               as.numeric(2*var(y)),
               10000)
  
  tauStore = array(NA, dim=c(length(sigmaBetaVec), length(sigmaVec), p))
  
  for (sb in 1 : length(sigmaBetaVec)) {
    for (sv in 1 : length(sigmaVec)) {
      for (j in 1 : p) {
        design = designMain2[,1:df2[[j]],j]
        
        PriorSigma = sigmaBetaVec[sb]*diag(df2[[j]])
        PriorSigmaInv <- chol2inv(chol(PriorSigma))
        mu = rep(0, df2[[j]])
        
        ## probability of being in top group
        muVar = solve(crossprod(design)  / sigmaVec[sv] +
                        PriorSigmaInv)
        
        logRATIOvec = rep(NA, 100)
        for (rr in 1 : length(logRATIOvec)) {
          yRandom = rnorm(n, mean=0, sd=sqrt(sigmaVec[sv]))
          
          muBetaRandom = muVar %*% (t(design) %*% yRandom/sigmaVec[sv] +
                                      PriorSigmaInv %*% mu)
          
          logRATIOvec[rr] = dmvn0Ander(mean=mu, sigma=PriorSigma) -
            dmvn0Ander(mean=muBetaRandom, sigma=muVar)
        }
        
        PROBvec = rep(NA, length(tauVec))
        for (tt in 1 : length(tauVec)) {
          tempTau = tauVec[tt]
          tempLogTau = logTauVec[tt]
          
          tauTimesRatio = exp(tempLogTau + logRATIOvec)
          PROBvec[tt] = mean(tauTimesRatio / ((1 - tempTau) + tauTimesRatio))
        }
        w = which(PROBvec < alphaMain)
        if (length(w) >= 1) {
          tauStore[sb,sv,j] = tauVec[w[length(w)]]
        } else {
          tauStore[sb,sv,j] = tauVec[1]
        }
      }
    }
  }
  
  sigmaBetaIntVec = c(0.0000001, seq(0.0001, 0.001, length=5),
                      seq(0.002, 0.01, length=5),
                      seq(0.02, 0.1, length=5),
                      seq(0.2, 1, length=5),
                      3,5,10,100,10000)
  
  tauIntStore = array(NA, dim=c(length(sigmaBetaIntVec), length(sigmaVec), choose(p,2)))
  
  for (sb in 1 : length(sigmaBetaIntVec)) {
    for (sv in 1 : length(sigmaVec)) {
      for (j in 1 : choose(p,2)) {
        j1 = intList[j,1]
        j2 = intList[j,2]
        
        design = designInt2[,1:dfI2[[j]],j]
        
        PriorSigmaInt = sigmaBetaIntVec[sb]*diag(dfI2[[j]])
        PriorSigmaIntInv <- chol2inv(chol(PriorSigmaInt))
        muInt = rep(0, dfI2[[j]])
        
        ## probability of being in top group
        muVar = solve(crossprod(design)  / sigmaVec[sv] +
                        PriorSigmaIntInv)
        
        
        logRATIOvec = rep(NA, 100)
        for (rr in 1 : length(logRATIOvec)) {
          yRandom = rnorm(length(y), mean=0, sd=sqrt(sigmaVec[sv]))
          
          muBetaRandom = muVar %*% (t(design) %*% yRandom/sigmaVec[sv] +
                                      PriorSigmaIntInv %*% muInt)
          
          logRATIOvec[rr] = dmvn0Ander(mean=muInt, sigma=PriorSigmaInt) -
            dmvn0Ander(mean=muBetaRandom, sigma=muVar)
        }
        PROBvec = rep(NA, length(tauVec))
        for (tt in 1 : length(tauVec)) {
          tempTau = tauVec[tt]
          tempLogTau = logTauVec[tt]
          
          tauTimesRatio = exp(tempLogTau + logRATIOvec)
          
          PROBvec[tt] = mean(tauTimesRatio / ((1 - tempTau) + tauTimesRatio))
        }
        
        w = which(PROBvec < alphaInt)
        if (length(w) >= 1) {
          tauIntStore[sb,sv,j] = tauVec[w[length(w)]]
        } else {
          tauIntStore[sb,sv,j] = tauVec[1]
        }
      }
    }
  }
  
  likelihood = array(NA, dim=c(nChains, nScans, n))
  
  print("Beginning MCMC")
  
  for (ni in 2 : nScans) {
    for (nc in 1 : nChains) {
      if (nc == 1 & ni %% 100 == 0) print(paste(ni, "MCMC scans have finished"))
      
      ## update sigma
      aStar = aSig + n/2
      bStar = bSig + sum((y - (designC %*% betaCPost[nc,ni-1,]) - apply(TempLinPred[nc,,], 1, sum) -
                            apply(TempLinPredInt[nc,,], 1, sum))^2)/2
      
      sigmaPost[nc,ni] = 1/rgamma(1,aStar,bStar)
      
      sigmaBetaPost[nc,ni] = min(1/rgamma(1, aSigBeta + sum(unlist(df2)*gammaPost[nc,ni-1,])/2,
                                          bSigBeta + sum(betaPost[nc,ni-1,,]^2)/(2)), 500)
      
      ## update sigmaBetaInt
      sigmaBetaIntPost[nc,ni] = min(1/rgamma(1, aSigBeta + sum(unlist(dfI2)*gammaIntPost[nc,ni-1,])/2,
                                             bSigBeta + sum(betaIntPost[nc,ni-1,,]^2)/(2)), 500)
      
      ## Update gamma and beta
      
      ## Update regression coefficients and variable inclusion parameters
      tempBeta = betaPost[nc,ni-1,,]
      for (j in 1 : p) {
        yStar = y - (designC %*% betaCPost[nc,ni-1,]) - apply(TempLinPred[nc,,-j], 1, sum) -
          apply(TempLinPredInt[nc,,], 1, sum)
        
        design = designMain2[,1:df2[[j]],j]
        
        PriorSigma = sigmaBetaPost[nc,ni]*diag(df2[[j]])
        PriorSigmaInv <- chol2inv(chol(PriorSigma))
        mu = rep(0, df2[[j]])
        
        ## probability of being in top group
        muVar = chol2inv(chol(crossprod(design)  / sigmaPost[nc,ni] + PriorSigmaInv))
        
        ## Fix tau as a function of other parameters
        wSB = which(sigmaBetaVec < sigmaBetaPost[nc,ni])
        wSV = which(sigmaVec > sigmaPost[nc,ni])
        
        if (length(wSB) == 0 | length(wSV) == 0) {
          tau = 0
        } else {
          tau = tauStore[wSB[length(wSB)],wSV[1],j]
        }
        
        ## probability of being in group zero
        p0 = log(1 - tau)
        
        muBeta = muVar %*% (t(design) %*% yStar/sigmaPost[nc,ni] + PriorSigmaInv %*% mu)
        p1 = log(tau) + dmvn0Ander(mean=mu, sigma=PriorSigma) - dmvn0Ander(mean=muBeta, sigma=muVar)
        
        maxlog = max(p0,p1)
        
        p0new = exp(-maxlog + p0)
        p1new = exp(-maxlog + p1)
        
        gammaPost[nc,ni,j] = sample(0:1, size=1, prob=c(p0new,p1new))
        if (VariableSel == FALSE) {
          gammaPost[nc,ni,j] = 1
        }
        
        tempBeta[j,] = rep(0, dfMain)
        if (gammaPost[nc,ni,j] == 1) {
          tempBeta[j,1:df2[[j]]] = rmnvAnder(mean=muBeta, sigma=muVar)
        }
        
        TempLinPred[nc,,j] = design %*% tempBeta[j,1:df2[[j]]]
      }
      betaPost[nc,ni,,] = tempBeta
      
      ## Update regression coefficients and variable inclusion parameters interactions
      tempBetaInt = betaIntPost[nc,ni-1,,]
      for (j in 1 : choose(p,2)) {
        j1 = intList[j,1]
        j2 = intList[j,2]
        
        yStar = y - (designC %*% betaCPost[nc,ni-1,]) - apply(TempLinPred[nc,,], 1, sum) -
          apply(TempLinPredInt[nc,,-j], 1, sum)
        
        design = designInt2[,1:dfI2[[j]],j]
        
        PriorSigmaInt = sigmaBetaIntPost[nc,ni]*diag(dfI2[[j]])
        PriorSigmaIntInv <- chol2inv(chol(PriorSigmaInt))
        muInt = rep(0, dfI2[[j]])
        
        ## probability of being in top group
        muVar = chol2inv(chol(crossprod(design)  / sigmaPost[nc,ni] + PriorSigmaIntInv))
        
        ## Fix tau as a function of other parameters
        wSB = which(sigmaBetaIntVec < sigmaBetaIntPost[nc,ni])
        wSV = which(sigmaVec > sigmaPost[nc,ni])
        
        if (length(wSB) == 0 | length(wSV) == 0) {
          tau = 0
        } else {
          tau = tauIntStore[wSB[length(wSB)],wSV[1],j]
        }
        
        ## probability of being in group zero
        p0 = log(1 - tau)
        
        muBeta = muVar %*% (t(design) %*% yStar/sigmaPost[nc,ni] +
                              PriorSigmaIntInv %*% muInt)
        p1 = log(tau) + dmvn0Ander(mean=muInt, sigma=PriorSigmaInt) - dmvn0Ander(mean=muBeta, sigma=muVar)
        
        maxlog = max(p0,p1)
        
        p0new = exp(-maxlog + p0)
        p1new = exp(-maxlog + p1)
        
        gammaIntPost[nc,ni,j] = sample(0:1, size=1, prob=c(p0new,p1new))
        
        tempBetaInt[j,] = rep(0, dfInt)
        TempLinPredInt[nc,,j] = rep(0, n)
        if (gammaIntPost[nc,ni,j] == 1) {
          tempBetaInt[j,1:dfI2[[j]]] = rmnvAnder(mean=muBeta, sigma=muVar)
          
          TempLinPredInt[nc,,j] = design %*% tempBetaInt[j,1:dfI2[[j]]]
        }
      }
      betaIntPost[nc,ni,,] = tempBetaInt
      
      ## Update intercept and covariate parameters
      yStar = y - apply(TempLinPred[nc,,], 1, sum) -
        apply(TempLinPredInt[nc,,], 1, sum)
      
      muBetaC = chol2inv(chol((CC)/sigmaPost[nc,ni-1] + SigmaCInv)) %*%
        ((t(designC) %*% yStar)/sigmaPost[nc,ni-1])
      covBetaC = chol2inv(chol((CC)/sigmaPost[nc,ni-1] + SigmaCInv))
      
      betaCPost[nc,ni,] = rmnvAnder(mean=muBetaC, sigma = covBetaC)
      
      linPred = (designC %*% betaCPost[nc,ni,]) + apply(TempLinPred[nc,,], 1, sum) +
        apply(TempLinPredInt[nc,,], 1, sum)
      likelihood[nc,ni,] = (1/(sqrt(2*pi*sigmaPost[nc,ni]))) *
        exp((-(y - linPred)^2)/(2*sigmaPost[nc,ni]))
    }
  }
  
  for (nc in 1 : nChains) {
    for (ni in 1 : nScans) {
      
      for (j in 1 : p) {
        ## re-scale the coefficients back to their original scale
        betaPost[nc,ni,j,1:df2[[j]]] =  
          as.vector(QmatMain2[[j]] %*% diag(1/sqrt(DvecMain2[[j]])) %*% betaPost[nc,ni,j,1:df2[[j]]])
        
        betaPost[nc,ni,j,1] = betaPost[nc,ni,j,1] / sdVecMain[1,j]
        
        if (df2[[j]] == 2) {
          betaPost[nc,ni,j,2:(df[[j]] + 1)] = 
            as.vector(QmatMain[[j]] %*% (1/sqrt(DvecMain[[j]])) %*% 
                        betaPost[nc,ni,j,2:df2[[j]]]) / sdVecMain[2:(df[[j]] + 1),j]
        } else {
          betaPost[nc,ni,j,2:(df[[j]] + 1)] = 
            as.vector(QmatMain[[j]] %*% diag(1/sqrt(DvecMain[[j]])) %*% 
                        betaPost[nc,ni,j,2:df2[[j]]]) / sdVecMain[2:(df[[j]] + 1),j] 
        }
        
        betaCPost[nc,ni,1] = betaCPost[nc,ni,1] -
          as.numeric(meanVecMain[1:(df[[j]] + 1),j] %*% betaPost[nc,ni,j,1:(df[[j]] + 1)])
      }
      
      for (j in 1 : choose(p,2)) {
        betaIntPost[nc,ni,j,1:dfI2[[j]]] = 
          as.vector(QmatInt2[[j]] %*% diag(1/sqrt(DvecInt2[[j]])) %*% betaIntPost[nc,ni,j,1:dfI2[[j]]])
        
        betaIntPost[nc,ni,j,1] = betaIntPost[nc,ni,j,1] / sdVecInt[1,j]
        
        if (dfI2[[j]] == 2) {
          betaIntPost[nc,ni,j,2:dfI[[j]]] = 
            as.vector(QmatInt[[j]] %*% (1/sqrt(DvecInt[[j]])) %*% 
                        betaIntPost[nc,ni,j,2:dfI2[[j]]]) / sdVecInt[2:dfI[[j]],j]
        } else {
          betaIntPost[nc,ni,j,2:dfI[[j]]] = 
            as.vector(QmatInt[[j]] %*% diag(1/sqrt(DvecInt[[j]])) %*% 
                        betaIntPost[nc,ni,j,2:dfI2[[j]]]) / sdVecInt[2:dfI[[j]],j]
        }
        
        betaCPost[nc,ni,1] = betaCPost[nc,ni,1] -
          as.numeric(meanVecInt[1:dfI[[j]],j] %*% betaIntPost[nc,ni,j,1:dfI[[j]]])
        
      }
    }
  }
  
  keep = seq(nBurn + 1, nScans, by=thin)
  
  WAIC = -2*(sum(log(apply(likelihood[,keep,], 3, mean))) -
               sum(apply(log(likelihood[,keep,]), 3, sd)^2))
  
  gammaInt = data.frame(intList)
  gammaInt$PIP = apply(gammaIntPost[,keep,], 3, mean)
  row.names(gammaInt) = NULL

  l = list(gamma = apply(gammaPost[,keep,], 3, mean),
           gammaInt = gammaInt,
           betaPost = betaPost[,keep,,],
           betaIntPost = betaIntPost[,keep,,],
           BasisFunctions = BasisFunctions,
           BasisFunctionsInt = BasisFunctionsInt,
           t=t, df=df, dfI = dfI, p=p, WAIC=WAIC)

  return(l)
}

