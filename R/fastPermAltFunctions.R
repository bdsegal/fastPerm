
xiFun <- function(beta,m){
  beta[1]*m + beta[2]/(1+exp(-beta[3]*m))
}

xiFunLin <- function(beta,m){
    beta*m
}

optFunNorm <- function(beta, m, pHat){
  sum((log(pHat) - log(pnorm(xiFun(beta,m), lower.tail=FALSE)))^2)
}

optFunT <- function(beta, m, pHat, nx, ny){
  sum((log(pHat) - log(pt(xiFun(beta,m), df=nx+ny, lower.tail=FALSE)))^2)
}

optFunNormLin <- function(beta, m, pHat){
  sum((log(pHat) - log(pnorm(xiFunLin(beta,m), lower.tail=FALSE)))^2)
}

optFunTLin <- function(beta, m, pHat, nx, ny){
  sum((log(pHat) - log(pt(xiFunLin(beta,m), df=nx+ny, lower.tail=FALSE)))^2)
}

  
# alternatives to Poisson regression
# proposed algorithm --------------------------------------
fastPermAlt <- function(x, y, testStat = ratioMean, B=1000, adjusted=FALSE){
#' Fast approximation of small permutation p-values
#'
#' This function approximates the p-value (two-sided) for a two
#' sample permutation test by partitioning the permutation space.
#' This function is most useful for small p-values, e.g. p < 10^-6.
#' @param x First sample
#' @param y Second sample
#' @param testStat Test statistic, defaults to the ratio of the
#' means (ratioMean). The difference in means (diffMean) is also available.
#' @param B Number of Monte Carlo iterations within each partition. Defaults to 1,000.
#' @keywords fast permtution test two sample
#' @export
#' @examples
#' x <- rexp(100, 5)
#' y <- rexp(100, 2)
#' mStopRatioMean(x, y)
#' fastPermAlt(x, y)
#'
#' x <- rnorm(110, 0, 1)
#' y <- rnorm(110, 1, 1)
#' mStopDiffMean(x, y)
#' fastPermAlt(x, y, testStat = diffMean)

  if (attributes(testStat)$comparison == "ratio" & (min(x,y)<0)) {
    stop("Some sample values <=0. With ratio statistics, all values must be greater than 0.")
  }
  
  nx <- length(x)
  ny <- length(y)
  
  # need nx <= ny. swap input if not correctly ordered
  if (ny<nx){
    xtemp <- x
    x <- y
    y <- xtemp
    
    nxtemp <- nx
    nx <- ny
    ny <- nxtemp
  }
   
  pmf <- PiLgamma(nx, ny)
  
  mMax <- as.numeric(names(pmf)[which(pmf == max(pmf))])
  t0 <- testStat(x, y)
  
  countTemp <-  rep(0, nx)
  names(countTemp) <- 1:nx
  m <- 1
  sumGTzero <- TRUE
  
  while (sumGTzero & (m <= mMax[1])) {

    tb <- rep(NA, B)
    
    for (b in 1:B){
      starX <- sample(1:nx, size = m, replace = FALSE)
      starY <- sample(1:ny, size = m, replace = FALSE)

      xNew <- c(x[-starX], y[starY])
      yNew <- c(y[-starY], x[starX])
      
      tb[b] <-  testStat(xNew, yNew)
    }
    
    countTemp[m] <- sum(tb >= t0)
    
    sumGTzero <- (countTemp[m] > 0)

    m <- m+1

  }

  mStop <- m-1
    
  # setup partitions for prediction, to be symmetric about mMax
  if (length(mMax)==1 & mMax[1]==nx){
    mPred <- 1:mMax
  } else if (length(mMax)==1 & mMax[1]<nx){
    mPred <- c(1:mMax,((mMax-1):0)[1:(nx-mMax)])
  } else if (length(mMax)==2 & mMax[2]==nx){
    mPred <- c(1:mMax[1], mMax[1])
  } else {
    mPred <- c(1:mMax[1], (mMax[1]:0)[1:(nx - mMax[1])])
  }
  
  # data frame for regression
  pHat <- (countTemp[1:(mStop-1)] + 1*adjusted)/B
  
  betaNorm <- optim(c(0.01,1,.5), optFunNorm, m=1:length(pHat), pHat=pHat)$par
  betaT <- optim(c(0.01,1,.5), optFunT, m=1:length(pHat), pHat=pHat, nx=nx, ny=ny)$par
  
  betaNormLin <- optimize(optFunNormLin, interval=c(0,1), m=1:length(pHat), pHat=pHat)$minimum
  betaTLin <- optimize(optFunTLin, interval=c(0,1), m=m, pHat=1:length(pHat), nx=nx, ny=ny)$minimum
  
  # betaNorm <- optim(c(.1, .1, .1), optFunNorm, m=m, pHat=pHat, method="CG", hessian=TRUE)
  # betaNorm <- optim(c(2,2,2), optFunNorm, m=m, pHat=pHat, method="L-BFGS-B", hessian=TRUE)

  # new data for predictions in partitions m = mPred
  xiNorm <- xiFun(beta=betaNorm, m=mPred)
  xiT <- xiFun(beta=betaT, m=mPred)
  xiNormLin <- xiFunLin(beta=betaNormLin, m=mPred)
  xiTLin <- xiFunLin(beta=betaTLin, m=mPred)
  
  # if nx < ny, set p-value in 0 partition to 1 and 
  # estimate p-value in partition nx
  if (nx != ny) { 
    pNorm <- c(1, pnorm(xiNorm, lower.tail=FALSE)) %*% pmf
    pT <- c(1, pt(xiT, df=nx+ny, lower.tail=FALSE)) %*% pmf
    pNormLin <- c(1, pnorm(xiNormLin, lower.tail=FALSE)) %*% pmf
    pTLin <- c(1, pt(xiTLin, df=nx+ny, lower.tail=FALSE)) %*% pmf
    
  # if nx==ny, set both the 0 and nx partition to 1
  } else { 
    pNorm <- c(1, pnorm(xiNorm[-length(xiNorm)], lower.tail=FALSE), 1) %*% pmf
    pT <- c(1, pt(xiT[-length(xiT)], df=nx+ny, lower.tail=FALSE), 1) %*% pmf
    pNormLin <- c(1, pnorm(xiNormLin[-length(xiNormLin)], lower.tail=FALSE), 1) %*% pmf
    pTLin <- c(1, pt(xiTLin[-length(xiTLin)], df=nx+ny, lower.tail=FALSE), 1) %*% pmf
  }

  ret <- list(pNorm = pNorm,
              pT = pT,
              pNormLin=pNormLin,
              pTLin=pTLin,
              mStop = mStop,
              betaNorm=betaNorm,
              betaT=betaT,
              betaNormLin=betaNormLin,
              betaTLin=betaTLin,
              xiNorm=xiNorm,
              xiNormLin=xiNormLin,
              xiT=xiT,
              xiTLin=xiTLin,
              mPred=mPred,
              B = B,
              t0 = t0,
              comparison = attributes(testStat)$comparison,
              summary = attributes(testStat)$summary)
        
  class(ret) <- "fastPermAlt"
  
  return(ret)
}

print.fastPermAlt <- function(fp){
  #' Print function for fastPerm
  #'
  #' This function prints the results of fastPerm
  #' @param fp Output from the fastPerm function
  #' @keywords fastPerm print
  #' @export
  #' @examples
  #' x <- rexp(100, 5)
  #' y <- rexp(100, 2)
  #' print(fastPerm(x, y, testStat = ratioMean))
  
  result <- paste("    fastPerm Two Sample Test   \n\n",
    fp$comparison, " of ", fp$summary, "s\n", 
    prettyNum(fp$B, big.mark = ","), " iterations within partitions",
    "\n\nobserved statistic = ", signif(fp$t0,3),
    "\n\np-value \t Approximation", 
    "\n-------------------------------",
    "\n", signif(fp$pT,3), "\t t ",
    "\n", signif(fp$pNorm,3), "\t normal",
    "\n", signif(fp$pTLin,3), "\t t, linear",
    "\n", signif(fp$pNormLin,3), "\t normal, linear xi",
    "\n\nmStop = ", fp$mStop,
    sep = "")
  
  writeLines(result)
}