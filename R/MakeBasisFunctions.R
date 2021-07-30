#' CreatePCA
#'
#' @param x matrix of exposures
#' @param PCAthresholdBasis number between 0 and 1
#'
#' @importFrom refund fpca.face

CreatePCA = function(x, PCAthresholdBasis) {
   BasisFunctions = list()
   BasisFunctionsInt = list()
   df = list()

   n = dim(x[[1]])[1]
   t = dim(x[[1]])[2]
   p = length(x)

   ## Find the degrees of freedom of basis functions for each exposure
   for (j in 1 : p) {
     covx = cov(x[[j]])
     pca = fpca.face(Y=covx, knots=10, pve=PCAthresholdBasis)

     w = dim(pca$efunctions)[2]
     df[[j]] = w
   }

   ## Do the same thing for interactions
   intList = expand.grid(1:p, 1:p)
   w = which(intList[,1] < intList[,2])
   intList = intList[w,]

   dfI = list()
   for (i in 1 : choose(p,2)) {
     j1 = intList[i,1]
     j2 = intList[i,2]
     dfI[[i]] = (df[[j1]]+1)*(df[[j2]]+1)
   }

   ## Find the maximum number of main and interaction basis functions
   dfMain = max(unlist(df)) + 1
   dfInt = prod(sort(unlist(df)+1, decreasing=TRUE)[1:2])

   ## Save the basis functions for the main effect DLAGs
   for (j in 1 : p) {
     covx = cov(x[[j]])
     pca = fpca.face(Y=covx, knots=10, pve=PCAthresholdBasis)

     w = dim(pca$efunctions)[2]
     BasisFunctions[[j]] = matrix(0, t, dfMain)
     BasisFunctions[[j]][,1] = rep(1,t)
     BasisFunctions[[j]][,2:(w+1)] = pca$efunctions[,1:w]
   }

   ## Save the basis functions for interaction surfaces
   for (i in 1 : choose(p,2)) {
     j1 = intList[i,1]
     j2 = intList[i,2]

     BasisFunctionsInt[[i]] = array(0, dim=c(t, t, dfInt))

     counter = 1
     for (k1 in 1 : (df[[j1]] + 1)) {
       for (k2 in 1 : (df[[j2]] + 1)) {
         for (t1 in 1 : t) {
           for (t2 in 1 : t) {
             BasisFunctionsInt[[i]][t1,t2,counter] = BasisFunctions[[j1]][t1,k1] *
               BasisFunctions[[j2]][t2,k2]
           }
         }
         counter = counter + 1
       }
     }
   }

   ## Now create the design matrices to be used for modeling
   designMain = array(NA, dim=c(n, dfMain, p))
   for (j in 1 : p) {
     design = matrix(NA, n, df[[j]] + 1)
     design[,1] = apply(x[[j]], 1, sum)
     for (k in 2 : (df[[j]] + 1)) {
       design[,k] = apply(t(t(x[[j]])*BasisFunctions[[j]][,k]), 1, sum)
     }
     designMain[,1:(df[[j]] + 1),j] = design
   }

   designInt = array(0, dim=c(n, dfInt, choose(p, 2)))
   for (i in 1 : choose(p,2)) {
     j1 = intList[i,1]
     j2 = intList[i,2]
     for (k in 1 : dfI[[i]]) {
       for (t1 in 1 : t) {
         for (t2 in 1 : t) {
           designInt[,k,i] = designInt[,k,i] + x[[j1]][,t1]*x[[j2]][,t2]*
             BasisFunctionsInt[[i]][t1,t2,k]
         }
       }
     }
   }


   l = list(df=df, dfI=dfI, dfMain=dfMain, dfInt=dfInt,
            designMain=designMain, designInt=designInt,
            BasisFunctions=BasisFunctions, BasisFunctionsInt=BasisFunctionsInt)

   return(l)
}

#' CreateCubic
#'
#' @param x matrix of exposures


CreateCubic = function(x) {
  BasisFunctions = list()
  BasisFunctionsInt = list()
  df = list()

  n = dim(x[[1]])[1]
  t = dim(x[[1]])[2]
  p = length(x)

  ## Find the degrees of freedom of basis functions for each exposure
  for (j in 1 : p) {
    df[[j]] = 3
  }

  ## Do the same thing for interactions
  intList = expand.grid(1:p, 1:p)
  w = which(intList[,1] < intList[,2])
  intList = intList[w,]

  dfI = list()
  for (i in 1 : choose(p,2)) {
    j1 = intList[i,1]
    j2 = intList[i,2]
    dfI[[i]] = (df[[j1]]+1)*(df[[j2]]+1)
  }

  ## Find the maximum number of main and interaction basis functions
  dfMain = max(unlist(df)) + 1
  dfInt = prod(sort(unlist(df)+1, decreasing=TRUE)[1:2])

  ## Save the basis functions for the main effect DLAGs
  for (j in 1 : p) {
    BasisFunctions[[j]] = matrix(0, t, dfMain)
    BasisFunctions[[j]][,1] = rep(1,t)
    BasisFunctions[[j]][,2] = (1:t)/t
    BasisFunctions[[j]][,3] = ((1:t)/t)^2
    BasisFunctions[[j]][,4] = ((1:t)/t)^3
  }

  ## Save the basis functions for interaction surfaces
  for (i in 1 : choose(p,2)) {
    j1 = intList[i,1]
    j2 = intList[i,2]

    BasisFunctionsInt[[i]] = array(0, dim=c(t, t, dfInt))

    counter = 1
    for (k1 in 1 : (df[[j1]] + 1)) {
      for (k2 in 1 : (df[[j2]] + 1)) {
        for (t1 in 1 : t) {
          for (t2 in 1 : t) {
            BasisFunctionsInt[[i]][t1,t2,counter] = BasisFunctions[[j1]][t1,k1] *
              BasisFunctions[[j2]][t2,k2]
          }
        }
        counter = counter + 1
      }
    }
  }

  ## Now create the design matrices to be used for modeling
  designMain = array(NA, dim=c(n, dfMain, p))
  for (j in 1 : p) {
    design = matrix(NA, n, df[[j]] + 1)
    design[,1] = apply(x[[j]], 1, sum)
    for (k in 2 : (df[[j]] + 1)) {
      design[,k] = apply(t(t(x[[j]])*BasisFunctions[[j]][,k]), 1, sum)
    }
    designMain[,1:(df[[j]] + 1),j] = design
  }

  designInt = array(0, dim=c(n, dfInt, choose(p, 2)))
  for (i in 1 : choose(p,2)) {
    j1 = intList[i,1]
    j2 = intList[i,2]
    for (k in 1 : dfI[[i]]) {
      for (t1 in 1 : t) {
        for (t2 in 1 : t) {
          designInt[,k,i] = designInt[,k,i] + x[[j1]][,t1]*x[[j2]][,t2]*
            BasisFunctionsInt[[i]][t1,t2,k]
        }
      }
    }
  }


  l = list(df=df, dfI=dfI, dfMain=dfMain, dfInt=dfInt,
           designMain=designMain, designInt=designInt,
           BasisFunctions=BasisFunctions, BasisFunctionsInt=BasisFunctionsInt)

  return(l)
}

#' CreateTPS
#'
#' @param x matrix of exposures
#' @importFrom npreg basis_tps

CreateTPS = function(x) {
  BasisFunctions = list()
  BasisFunctionsInt = list()
  df = list()

  n = dim(x[[1]])[1]
  t = dim(x[[1]])[2]
  p = length(x)

  nKnots = min(ceiling(t/2), 10)
  knots = seq(from=2, to= t-1, length=nKnots)

  knotsInt = as.matrix(expand.grid(knots, knots))

  ## Find the degrees of freedom of basis functions for each exposure
  for (j in 1 : p) {
    ## Add one for the linear term
    df[[j]] = nKnots + 1
  }

  ## Do the same thing for interactions
  intList = expand.grid(1:p, 1:p)
  w = which(intList[,1] < intList[,2])
  intList = intList[w,]

  dfI = list()
  for (i in 1 : choose(p,2)) {
    j1 = intList[i,1]
    j2 = intList[i,2]

    dataTPS = expand.grid(1:t, 1:t)
    Xtemp <- basis_tps(x=dataTPS, knots=knotsInt)

    dfI[[i]] = dim(Xtemp)[2] + 1
  }

  ## Find the maximum number of main and interaction basis functions
  dfMain = max(unlist(df)) + 1
  dfInt = max(unlist(dfI))

  ## Save the basis functions for the main effect DLAGs
  for (j in 1 : p) {
    dataTPS = expand.grid(1:t)

    Xtemp <- basis_tps(x=dataTPS, knots=knots)
    BasisFunctions[[j]] = matrix(0, t, dfMain)
    BasisFunctions[[j]][,1] = rep(1,t)
    BasisFunctions[[j]][,2 : dfMain] = Xtemp
  }

  ## Save the basis functions for interaction surfaces
  for (i in 1 : choose(p,2)) {
    j1 = intList[i,1]
    j2 = intList[i,2]

    BasisFunctionsInt[[i]] = array(0, dim=c(t, t, dfInt))

    dataTPS = expand.grid(1:t, 1:t)
    Xtemp <- basis_tps(x=dataTPS, knots=knotsInt)

    for (t1 in 1 : t) {
      for (t2 in 1 : t) {
        BasisFunctionsInt[[i]][t1,t2,1] = 1
        wJ = which(dataTPS$Var1 == t1 & dataTPS$Var2 == t2)

        BasisFunctionsInt[[i]][t1,t2,2:dfInt] = Xtemp[wJ,]

      }
    }
  }

  ## Now create the design matrices to be used for modeling
  designMain = array(NA, dim=c(n, dfMain, p))
  for (j in 1 : p) {
    design = matrix(NA, n, df[[j]] + 1)
    design[,1] = apply(x[[j]], 1, sum)
    for (k in 2 : (df[[j]] + 1)) {
      design[,k] = apply(t(t(x[[j]])*BasisFunctions[[j]][,k]), 1, sum)
    }
    designMain[,1:(df[[j]] + 1),j] = design
  }

  designInt = array(0, dim=c(n, dfInt, choose(p, 2)))
  for (i in 1 : choose(p,2)) {
    j1 = intList[i,1]
    j2 = intList[i,2]
    for (k in 1 : dfI[[i]]) {
      for (t1 in 1 : t) {
        for (t2 in 1 : t) {
          designInt[,k,i] = designInt[,k,i] + x[[j1]][,t1]*x[[j2]][,t2]*
            BasisFunctionsInt[[i]][t1,t2,k]
        }
      }
    }
  }


  l = list(df=df, dfI=dfI, dfMain=dfMain, dfInt=dfInt,
           designMain=designMain, designInt=designInt,
           BasisFunctions=BasisFunctions, BasisFunctionsInt=BasisFunctionsInt)

  return(l)
}
