# Implementation of SAMC algorithm 

swapFun <- function(x, y, L){
  #' Utility function for SAMC function
  #'
  #' This function is used by SAMC to produce the next random permutation
  #' @param x First sample
  #' @param y Second sample
  #' @param L Number of elements to swap
  #' @keywords SAMC utility
  #' @export
  #' @examples
  #' swapFun(x, y, L)

  nx <- length(x)
  ny <- length(y)
  
  starX <- sample(1:nx, L)
  starY <- sample(1:ny, L)
  
  xstar <- c(x[-starX], y[starY])
  ystar <- c(y[-starY], x[starX])
  
  return(list(xstar = xstar, ystar = ystar))
}

gammaFun <- function(b, b0){
  #' Utility function for SAMC function
  #'
  #' This function is used by SAMC to deteroriate the effect of updates.
  #' @param b Current iteration of the SAMC algorithim
  #' @param bo Iteration at which the updates begin to deteriorate.
  #' @keywords SAMC utility
  #' @export
  #' @examples
  #' gammaFun(b, b0)
  
  b0 / max(b, b0)
}

SAMC <- function(x, y, testStat=ratioMean, B=10e4, m=300, b0=1000){
#' Stochastic Approximation Monte Carlo (SAMC)
#'
#' This function implements the SAMC algorithm developed by Liang et al. (2007) and 
#' tailored to p-value estimation by Yu et al. (2011)
#' @param x First sample
#' @param y Second sample
#' @param testStat Test statistic, defaults to the ratio of the
#' means (ratioMean). Other choices are diffMean, ratioMedian, and
#' diffMedian.
#' @param B Total Number of Monte Carlo iterations. Defaults to 10e4.
#' @param m Total number of regions in the SAMC algorithm. Note: SAMC regions are 
#' different than fastPerm partitions.
#' @param b0 Iteration at which the updates begin to decay.
#' @keywords SAMC two sample
#' @export
#' @examples
#' x <- rexp(100, 5)
#' y <- rexp(100, 2)
#' sam <- SAMC(x, y)
#' sam
#'
#' x <- rnorm(110, 0, 1)
#' y <- rnorm(110, 1, 1)
#' sam <- SAMC(x, y, testStat = diffMean)
#' sam
#' plot(sam)

  if (attributes(testStat)$comparison == "ratio" & (min(x,y)<0)) {
    stop("Some sample values <=0. With ratio statistics, all values must be greater than 0.")
  }
  
  ### keep track of algorithm
  numAccepts <- 0
  piObs <- rep(0, m)
  
  ###########################
  # initialize quantitiess
  n <- min(length(x), length(y))
  
  t0 <- testStat(x, y)

  if (attributes(testStat)$comparison == "difference"){
    lambda <- seq(0, t0, length.out=m)
  } else {
    lambda <- seq(1, t0, length.out=m)
  }
  
  theta <- rep(0, m)
  piVec <- rep(1 / m, m)

  L <- ceiling(n * 0.075)
  
  # get initial test stat and partition
  tInit <- testStat(x, y)
  jInit <- sum(tInit >= lambda)
  
  for (b in 1:B){
    # propose update
    swapOut <- swapFun(x, y, L)
    xNew <- swapOut$xstar
    yNew <- swapOut$ystar

    # get new test stat and partition
    # T_new <- abs(t.test(xNew,yNew)$statistic)
    tNew <- testStat(xNew, yNew)
    jNew <- sum(tNew >= lambda)
    
    # get acceptance probability and decide whether to accept
    r = min(1, exp(theta[jInit] - theta[jNew]))
    accept <- (runif(n = 1) <= r)
      
    # set new (x,y) and update weights
    if (accept){
      x <- xNew
      y <- yNew
      
      I <- rep(0, m)
      I[jNew] <- 1
      theta <- theta + gammaFun(b, b0) * (I-piVec)
      
      tInit <- tNew
      jInit <- jNew
      
      numAccepts <- numAccepts + 1
      
    } else{
      I <- rep(0, m)
      I[jInit] <- 1
      theta <- theta + gammaFun(b, b0) * (I-piVec)
    }

    piObs <- piObs + I
  }

  # estimate p-value
  pval = exp(theta[m]) * piVec[m] / sum(exp(theta) * piVec)
  
  ret <- list(pval = pval,
    numAccepts = numAccepts,
    piObs = piObs,
    theta = theta,
    t0 = t0,
    B = B,
    comparison = attributes(testStat)$comparison,
    summary = attributes(testStat)$summary)

  class(ret) <- "SAMC"
  
  return(ret)
}

print.SAMC <- function(sam){
  #' Print function for SAMC
  #'
  #' This function prints the reslts and convergence diagnostics for the SAMC algorithm.
  #' If the SAMC converged, then it should have sampled nearly uniformly from all regions.
  #' The maximum discrepancy is the maximum difference between the number of times SAMC
  #' sampled from a region and the expected amount if the sampling were uniform.
  #' Yu et al. suggest that the maximum discrepancy should be less than 0.2.
  #' @param sam Output from the SAMC function
  #' @keywords SAMC summary convergence diagnostics
  #' @export
  #' @examples
  #' x <- rexp(100, 5)
  #' y <- rexp(100, 2)
  #' sam <- SAMC(x, y)
  #' print(sam)
  
  x <- sam$piObs
  m0 <- sum(x==0)
  p <- x/sum(x)
  m <- length(p)

  error <- (p-1/(m-m0))/(1/(m-m0))
  error <- error*(x>0)
  
  result <- paste("    SAMC Two Sample Test    \n\n",
    sam$comparison, " of ", sam$summary, "s\n",
    prettyNum(sam$B, big.mark=","), " total iterations",
    "\n\nobserved statistic = ", signif(sam$t0,3),
    "\np-value = ", signif(sam$pval,3),
    "\n\nmax discrepancy = ", signif(max(abs(error)),2),
    "\nalgorithm has converged if max discrepancy < 0.2", sep = "")
  
  writeLines(result)
}

plot.SAMC <- function(sam){
  #' plot function for SAMC
  #'
  #' This function plots the convergence diagnostics for the SAMC algorithm and prints
  #' out the associated values. If the SAMC converged, then it should have sampled
  #' nearly uniformly from all regions. The maximum discrepancy is the maximum difference
  #' between the number of times SAMC sampled from a region and the expected amount if the
  #' sampling were uniform. Yu et al. suggest that the maximum discrepancy should be less 
  #' than 0.2.
  #' @param sam Output from the SAMC function
  #' @keywords SAMC summary convergence diagnostics
  #' @export
  #' @examples
  #' x <- rexp(100, 5)
  #' y <- rexp(100, 2)
  #' sam <- SAMC(x, y)
  #' plot(sam)
  
  x <- sam$piObs
  m0 <- sum(x==0)
  p <- x/sum(x)
  m <- length(p)

  error <- (p-1/(m-m0))/(1/(m-m0))
  error <- error*(x>0)
  
  maxDiscrep <- signif(max(abs(error)),2)
  
  barplot(x, xlab="SAMC region", ylab="count", 
    main=paste("SAMC sampling distribution (uniform at convergence)",
    "\nmax discrepancy = ", maxDiscrep, " (<0.2 at convergence)",
    sep = ""))
}