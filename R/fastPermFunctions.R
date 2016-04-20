# pmf of partitions, and
# proposed algorithm

# pmf of partitions ---------------------------------------

Pi <- function(nx, ny){
  #' Exact pmf of partitions
  #'
  #' This function calculates the exact probability mass function
  #' of the partitions.
  #' Typically, this function is called from fastPerm.
  #' @param nx Size of first sample
  #' @param ny Size of second sample
  #' @keywords partition pmf
  #' @export
  #' @examples
  #' Pi(nx, ny)
  
  out <- rep(NA, min(nx, ny) + 1)
  names(out) <- 0:min(nx, ny)
  
  for (m in 0:min(nx, ny)){
    out[m + 1] <- choose(nx, m) * choose(ny, m) / choose(nx + ny, nx)
  }
  
  return(out)
}

PiLgammaNxSmaller <- function(nx, ny){
  #' Exact pmf of partitions with the log gamma function
  #' Assuming nx <= ny
  #'
  #' This function calculates the exact probability mass function
  #' of the partitions using the lgamma function.
  #' Typically, this function is called from fastPerm.
  #' @param nx Size of first sample
  #' @param ny Size of second sample
  #' @keywords partition pmf
  #' @export
  #' @examples
  #' PiLgamma(nx, ny)
  
  out <- rep(NA, min(nx, ny) + 1)
  names(out) <- 0:min(nx, ny)
  
  for (m in 0:min(nx, ny)){
    out[m + 1] <- exp(2*(lgamma(nx + 1) + lgamma(ny + 1) - lgamma(m + 1)) - 
      lgamma(nx + ny + 1) - lgamma(nx - m + 1) - lgamma(ny - m + 1))
  }
  
  return(out)
}

PiLgamma <- function(nx, ny){
  #' Exact pmf of partitions with the log gamma function
  #'
  #' This function calculates the exact probability mass function
  #' of the partitions using the lgamma function.
  #' Typically, this function is called from fastPerm.
  #' @param nx Size of first sample
  #' @param ny Size of second sample
  #' @keywords partition pmf
  #' @export
  #' @examples
  #' PiLgamma(nx, ny)
  
  nMin <- min(nx, ny)
  out <- rep(NA, nMin + 1)
  names(out) <- 0:min(nx, ny)
  
  for (m in 0:nMin){
    out[m + 1] <- exp(
    	lgamma(nx + 1) + lgamma(ny + 1) + lgamma(nx + ny - nMin + 1) +
    	lgamma(nMin + 1) - lgamma(nx - m + 1) - lgamma(ny - m + 1) -
    	lgamma(nx + ny +1) - 2*lgamma(m+1)
    	)
  }
  
  return(out)
}
# test statistics, used as input to fastPerm ---------------

diffMean <- function(x, y){
  #' Utility function for fastPerm and SAMC
  #'
  #' Difference in means test statistic used by fastPerm
  #' @param x First sample 
  #' @param y Second sample
  #' @keywords test statistic
  #' @export
  #' @examples
  #' diffMean(x, y )
  
  abs(mean(x) - mean(y))
}
attributes(diffMean) <- list(summary="mean", comparison="difference")

diffMedian <- function(x, y){
  #' Utility function for fastPerm and SAMC
  #'
  #' Difference in medians test statistic used by fastPerm.
  #' Note: In this version of fastPerm, median statistic may not be reliable.
  #' @param x First sample 
  #' @param y Second sample
  #' @keywords test statistic
  #' @export
  #' @examples
  #' diffMedian(x, y )
  
  abs(median(x) - median(y))
}
attributes(diffMedian) <- list(summary="median", comparison="difference")

ratioMean <- function(x, y){
  #' Utility function for fastPerm and SAMC
  #'
  #' Ratio of means test statistic used by fastPerm
  #' @param x First sample 
  #' @param y Second sample
  #' @keywords test statistic
  #' @export
  #' @examples
  #' ratioMean(x, y )
  
  max(mean(x) / mean(y), mean(y) / mean(x))
}
attributes(ratioMean) <- list(summary="mean", comparison="ratio")

ratioMedian <- function(x, y){
  #' Utility function for fastPerm and SAMC
  #'
  #' Ratio of medians test statistic used by fastPerm.
  #' Note: In this version of fastPerm, median statistic may not be reliable.
  #' @param x First sample 
  #' @param y Second sample
  #' @keywords test statistic
  #' @export
  #' @examples
  #' ratioMedian(x, y )
  
  max(median(x) / median(y), median(y) / median(x))
}
attributes(ratioMedian) <- list(summary="median", comparison="ratio")

# proposed algorithm --------------------------------------
fastPerm <- function(x, y, testStat = ratioMean, B=1000, adjusted=FALSE){
#' Fast approximation of small permutation p-values
#'
#' This function approximates the p-value (two-sided) for a two
#' sample permutation test by partitioning the permutation space.
#' This function is most useful for small p-values, e.g. p < 10^-6.
#' @param x First sample
#' @param y Second sample
#' @param testStat Test statistic, defaults to the ratio of the
#' means (ratioMean). Other choices are diffMean, ratioMedian, and
#' diffMedian. In the current version, the median statistics are experimental and may
#' not be reliable.
#' @param B Number of Monte Carlo iterations within each partition. Defaults to 1,000.
#' @keywords fast permtution test two sample
#' @export
#' @examples
#' x <- rexp(100, 5)
#' y <- rexp(100, 2)
#' mStopRatioMean(x, y)
#' fastPerm(x, y)
#'
#' x <- rnorm(110, 0, 1)
#' y <- rnorm(110, 1, 1)
#' mStopDiffMean(x, y)
#' fastPerm(x, y, testStat = diffMean)

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
  
  # while (sumGTzero & (m <= mMax[1])) {
  while (sumGTzero & (m <= nx)) {

    tb <- rep(NA, B)
    
    for (b in 1:B){
      starX <- sample(1:nx, size = m, replace = FALSE)
      starY <- sample(1:ny, size = m, replace = FALSE)

      xNew <- c(x[-starX], y[starY])
      yNew <- c(y[-starY], x[starX])
      
      tb[b] <-  testStat(xNew, yNew)
    }
    
    countTemp[m] <- sum(tb >= t0)
    
    if (countTemp[m] == 0){
      sumGTzero <- FALSE
    } else{
      m <- m+1
    }
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
  count <- c(B, countTemp[1:mStop]) + 1*adjusted
  countForReg <- data.frame(count = count, mReg = 0:mStop)
  
  poisFit <- glm(count ~ mReg, family = poisson, data = countForReg)
  
  # new data for predictions in partitions m = mPred
  mNewData <- data.frame(mReg = mPred)
  
  # if nx < ny, set p-value in 0 partition to 1 and 
  # estimate p-value in partition nx
  if (nx < ny) { 
      
    # note: predictions with type="response" are not reliable for small values. 
    # There is not a problem for type="link"
    pPoisCount <- c(B + 1*adjusted,
      exp(predict(poisFit, newdata=mNewData, type="link")))

    pwTilde <- pPoisCount %*% pmf / (B + 1*adjusted)
  
  # if nx==ny, set both the 0 and nx partition to 1
  } else { 
  
    mNewData <-data.frame(mReg = mNewData$mReg[-nx])
    
    pPoisCount <- c(B + 1*adjusted,
      exp(predict(poisFit, newdata = mNewData, type="link")),
      B + 1*adjusted)
    
    pwTilde <- pPoisCount %*% pmf / (B + 1*adjusted)

  }
  
  glmSummary <- summary(poisFit)

  ret <- list(pwTilde = pwTilde,
    mStop = mStop,
    deviance = glmSummary$deviance,
    aic = glmSummary$aic,
    df.residual = glmSummary$df.residual,
    B = B,
    t0 = t0,
    comparison = attributes(testStat)$comparison,
    summary = attributes(testStat)$summary)
        
  class(ret) <- "fastPerm"
  
  return(ret)
}

print.fastPerm <- function(fp){
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
    "\np-value = ", signif(fp$pwTilde,3),
    "\nmStop = ", fp$mStop, ", deviance = ", signif(fp$deviance,3), ", AIC = ",
    signif(fp$aic,3), sep = "")
  
  writeLines(result)
}