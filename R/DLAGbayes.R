#' Estimate effects of repeated exposures using distributed lag models
#'
#'
#' @param y              The outcome to be analyzed
#' @param x              A list where each element is a matrix containing the values
#'                       of each exposure over time
#' @param c              An n by q matrix of additional, time invariant covariates to adjust for
#' @param nScans         The number of MCMC scans to run
#' @param nBurn          The number of MCMC scans that will be dropped as a burn-in
#' @param thin           This number represents how many iterations between each scan
#'                       that is kept
#' @param PCAthresholdBasis   This number represents what percentage of variation we want our basis functions
#'                       from PCA to capture from the original exposures
#' @param PCAthresholdMain ?
#' @param PCAthresholdInt ?
#' @param VariableSel    Whether to perform variable selection on main effect distributed lag surfaces
#' @param Basis          Indicates what type of basis functions to use. The default is PCA, which then
#'                       uses the PCA threshold variable to determine how many basis functions to keep.
#'                       The user may also input "Cubic" to get Cubic polynomial basis functions
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom stats lm rgamma runif rbeta cov
#' @return A list of values that contain the posterior inclusion probabilities for the main effect
#'         distributed lags and interaction distributed lag surfaces. The list also contains the
#'         posterior distribution for the distributed lag parameters as well as the basis functions
#'         used to estimate their distributed lag function. These posterior distributions do not
#'         represent the distributed lag surfaces themselves, but these can be obtained using the
#'         other functions in the package.
#'
#' @export


DLAGbayes = function(y, x, c=NULL,
                     nScans = 2000, nBurn = 1000, thin = 2,
                     PCAthresholdBasis = 0.95,
                     PCAthresholdMain = 1,
                     PCAthresholdInt = 0.999,
                     VariableSel = TRUE, Basis = "PCA") {

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

  SigmaC = diag(10000, q+1)

  intList = expand.grid(1:p, 1:p)
  w = which(intList[,1] < intList[,2])
  intList = intList[w,]

  ## Create basis functions for main effects and interactions
  if (Basis == "PCA") {
    createMats = CreatePCA(x=x, PCAthresholdBasis = PCAthresholdBasis)
  } else {
    createMats = CreateCubic(x=x)
  }

  designMain = createMats$designMain
  designInt = createMats$designInt
  BasisFunctions = createMats$BasisFunctions
  BasisFunctionsInt = createMats$BasisFunctionsInt
  df = createMats$df
  dfI = createMats$dfI
  dfMain = createMats$dfMain
  dfInt = createMats$dfInt

  ## Scale the data appropriately
  meanVecMain = apply(designMain, 2:3, mean, na.rm=TRUE)
  meanVecInt = apply(designInt, 2:3, mean, na.rm=TRUE)

  for (j in 1 : p) {
    for (k in 1 : (df[[j]] + 1)) {
      designMain[,k,j] = (designMain[,k,j] - meanVecMain[k,j])
    }
  }

  for (j in 1 : choose(p, 2)) {
    for (k in 1 : dfI[[j]]) {
      designInt[,k,j] = (designInt[,k,j] - meanVecInt[k,j])
    }
  }

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
      ggModel = cv.glmnet(y=y, x=overallDesign)
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
      ggModel = cv.glmnet(y=y, x=cbind(overallDesign, as.matrix(c)))
      sigmaSqEst = mean((y - cbind(rep(1, n), overallDesign, as.matrix(c)) %*%
                           as.numeric(coef(ggModel, s="lambda.min")))^2)

      y = y/sqrt(sigmaSqEst)
      designMain = designMain/sqrt(sigmaSqEst)
      designInt = designInt/sqrt(sigmaSqEst)
      c = as.matrix(c)/sqrt(sigmaSqEst)
    }
  }

  QmatMain = DvecMain = df2 = list()
  designMain2 = array(NA, dim=dim(designMain))
  for (j in 1 : p) {
    tempX = designMain[,1 : (df[[j]] + 1),j]
    SVD = svd((1/n)*t(tempX) %*% tempX)
    nComponents = which((cumsum(SVD$d) / sum(SVD$d)) >= PCAthresholdMain)[1]
    df2[[j]] = nComponents
    QmatMain[[j]] = SVD$u[,1:nComponents]
    DvecMain[[j]] = SVD$d[1:nComponents]

    designMain2[,1:nComponents,j] = tempX %*% QmatMain[[j]] %*% diag(1/sqrt(DvecMain[[j]]))
  }

  QmatInt = DvecInt = dfI2 = list()
  designInt2 = array(NA, dim=dim(designInt))
  for (j in 1 : choose(p,2)) {
    tempX = designInt[,1:dfI[[j]],j]
    SVD = svd((1/n)*t(tempX) %*% tempX)
    nComponents = which((cumsum(SVD$d) / sum(SVD$d)) >= PCAthresholdInt)[1]
    dfI2[[j]] = nComponents
    QmatInt[[j]] = SVD$u[,1:nComponents]
    DvecInt[[j]] = SVD$d[1:nComponents]

    designInt2[,1:nComponents,j] =
      tempX %*% QmatInt[[j]] %*% diag(1/sqrt(DvecInt[[j]]))
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
  tauPost = tauIntPost = array(NA, dim=c(nChains, nScans))

  ## TODO MAKE THESE RANDOM
  gammaPost[,1,] = 0
  gammaIntPost[,1,] = 0
  betaPost[,1,,] = 0
  betaCPost[,1,] = 0
  betaIntPost[,1,,] = 0
  sigmaPost[,1] = rgamma(nChains, 1, 1)
  sigmaBetaPost[,1] = rgamma(nChains, 1, 1)
  sigmaBetaIntPost[,1] = rgamma(nChains, 1, 1)
  tauPost[,1] = runif(2)
  tauIntPost[,1] = runif(2)

  aSig = 0.001
  bSig = 0.001

  aTau = 3
  bTau = p

  aSigBeta = 0.5
  bSigBeta = 0.5

  TempLinPred = array(0, dim=c(nChains, n, p))
  TempLinPredInt = array(0, dim=c(nChains, n, choose(p,2)))

  likelihood = array(NA, dim=c(nChains, nScans, n))

  print("Starting MCMC")

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

      ## Update tau
      tauPost[nc,ni] = rbeta(1, aTau + sum(gammaPost[nc,ni-1,] == 1),
                             bTau + sum(gammaPost[nc,ni-1,] == 0))

      ## update tau int
      tauIntPost[nc,ni] = rbeta(1, aTau + sum(gammaIntPost[nc,ni-1,] == 1),
                                bTau + sum(gammaIntPost[nc,ni-1,] == 0))

      ## Update gamma and beta

      ## Update regression coefficients and variable inclusion parameters
      tempBeta = betaPost[nc,ni-1,,]
      for (j in 1 : p) {
        yStar = y - (designC %*% betaCPost[nc,ni-1,]) - apply(TempLinPred[nc,,-j], 1, sum) -
          apply(TempLinPredInt[nc,,], 1, sum)

        ## probability of being in group zero
        p0 = log(1 - tauPost[nc,ni])

        design = designMain2[,1:df2[[j]],j]

        PriorSigma = sigmaBetaPost[nc,ni]*diag(df2[[j]])
        mu = rep(0, df2[[j]])

        ## probability of being in top group
        muVar = solve(t(design) %*% design  / sigmaPost[nc,ni] +
                        solve(PriorSigma))
        muBeta = muVar %*% (t(design) %*% yStar/sigmaPost[nc,ni] +
                              solve(PriorSigma) %*% mu)
        p1 = log(tauPost[nc,ni]) + dmvnorm(rep(0, df2[[j]]), mean=mu, sigma=PriorSigma, log=TRUE) -
          dmvnorm(rep(0, df2[[j]]), mean=muBeta, sigma=muVar, log=TRUE)

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

        ## probability of being in group zero
        p0 = log(1 - tauIntPost[nc,ni])

        design = designInt2[,1:dfI2[[j]],j]

        PriorSigmaInt = sigmaBetaIntPost[nc,ni]*diag(dfI2[[j]])
        muInt = rep(0, dfI2[[j]])

        ## probability of being in top group
        muVar = solve(t(design) %*% design  / sigmaPost[nc,ni] +
                        solve(PriorSigmaInt))
        muBeta = muVar %*% (t(design) %*% yStar/sigmaPost[nc,ni] +
                              solve(PriorSigmaInt) %*% muInt)
        p1 = log(tauIntPost[nc,ni]) +
          dmvnorm(rep(0, dfI2[[j]]), mean=muInt,
                           sigma=PriorSigmaInt, log=TRUE) -
          dmvnorm(rep(0, dfI2[[j]]), mean=muBeta,
                           sigma=muVar, log=TRUE)

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

      muBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC)) %*%
        ((t(designC) %*% yStar)/sigmaPost[nc,ni-1])
      covBetaC = solve((t(designC) %*% designC)/sigmaPost[nc,ni-1] + solve(SigmaC))

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
        betaCPost[nc,ni,1] = betaCPost[nc,ni,1] -
          meanVecMain[1:(df[[j]] + 1),j] %*% QmatMain[[j]] %*%
          diag(1/sqrt(DvecMain[[j]])) %*% betaPost[nc,ni,j,1:df2[[j]]]
        betaPost[nc,ni,j,1:(df[[j]] + 1)] = QmatMain[[j]] %*% diag(1/sqrt(DvecMain[[j]])) %*%
          betaPost[nc,ni,j,1:df2[[j]]]
      }

      for (j in 1 : choose(p,2)) {
        betaCPost[nc,ni,1] = betaCPost[nc,ni,1] -
          meanVecInt[1:dfI[[j]],j] %*% QmatInt[[j]] %*% diag(1/sqrt(DvecInt[[j]])) %*% betaIntPost[nc,ni,j,1:dfI2[[j]]]
        betaIntPost[nc,ni,j,1:dfI[[j]]] = QmatInt[[j]] %*% diag(1/sqrt(DvecInt[[j]])) %*% betaIntPost[nc,ni,j,1:dfI2[[j]]]
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




