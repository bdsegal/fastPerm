# Functions derived in the Appendix -----------------------

muw <- function(m, xbar, ybar) {
  #' Utitlity function for mStop
  #'
  #' Calculates expected value of W in partition m
  #' @export
  
  m * (ybar - xbar)
}

g <- function(m, nx, ny, xbar, ybar){
  #' Utitlity function for mStop
  #'
  #' Calculates g function of the expected value of W
  #' @export
  
  numerator <- nx * xbar + muw(m = m, xbar = xbar, ybar = ybar)
  denominator <- ny * ybar - muw(m = m, xbar = xbar, ybar = ybar)
  return(ny / nx * numerator / denominator)
}

gp <- function(m, nx, ny, xbar, ybar){
  #' Utitlity function for mStop
  #'
  #' Calculates derivative of the g function of the expected value of W
  #' @export
  
  numerator <- ny * ybar + nx * xbar
  denominator <- (ny * ybar - muw(m = m, xbar = xbar, ybar = ybar))^2
  return(ny / nx * numerator / denominator)
}

gLin <- function(m, nx, ny, xbar, ybar){
  #' Utitlity function for mStop
  #'
  #' Calculates g function of the expected value of W for T = |x - y|
  #' @export
  return(xbar - ybar + (1 / nx + 1 / ny) * muw(m = m, xbar = xbar, ybar = ybar))
  
}

gpLin <- function(nx, ny){
  #' Utitlity function for mStop
  #'
  #' Calculates g' function of the expected value of W for T = |x - y|
  #' @export
  
  return((1 / nx + 1 / ny))
}

VarW <- function(m, nx, ny, x, y){
  #' Utitlity function for mStop
  #'
  #' Calculates variance of W
  #' @export
  
  varX <- var(x)
  varY <- var(y)
  return(m * (ny - m) / ny * varY + m * (nx - m) / nx * varX)
}

VarR <- function(m, nx, ny, x, y, xbar, ybar){
  #' Utitlity function for mStop
  #'
  #' Calculates variance of the ratio stat R
  #' @export
  
  return(gp(m, nx, ny, xbar, ybar)^2 * VarW(m, nx, ny, x, y))
}

VarRLin <- function(m,nx,ny, x, y, xbar, ybar){
  #' Utitlity function for mStop
  #'
  #' Calculates variance of the linear stat T = |x - y|
  #' @export
  
  return(gpLin(nx, ny)^2 * VarW(m, nx, ny, x, y))
}


xiFunRatioMean <- function(m, nx, ny, x, y){
  #' Utitlity function for mStop
  #'
  #' Calculates xi(m)
  #' @export
  
  N <- nx + ny
  # lambda <- nx / N
  
  xbar <- mean(x)
  ybar <- mean(y)
  
  Emean <- g(m, nx, ny, xbar, ybar)
  Evar <- VarR(m, nx, ny, x, y, xbar, ybar)
  
  t <- xbar / ybar

  xi <- (t - Emean) / sqrt(Evar)
  
  return(xi)
}

mStopRatioMean <- function(x, y, B=1000, plot = FALSE){
  #' Estimate partition at which fastPerm will stop for ratio statistics
  #'
  #' This function estimates the partition at which fastPerm will
  #' will stop, based on asymptotic approximation.
  #' @param x First sample
  #' @param y Second sample
  #' @param B Number of Monte Carlo iterations within each partition. Defaults to 1,000.
  #' @keywords fast permtution stopping partition
  #' @export
  #' @examples
  #' x <- rexp(100, 5)
  #' y <- rexp(100, 2)
  #' mStopRatioMean(x, y)

  # swap x and y if xbar < ybar
  if (mean(x) < mean(y)) {
    ytemp <- y
    y <- x
    x <- ytemp
  }

  nx <- length(x)
  ny <- length(y)
  
  pmf <- PiLgamma(nx, ny)
  mMax <- min(which(pmf == max(pmf)))

  u <- qnorm(1-1/B)
  m <- 1:min(mMax,(min(nx, ny)-1))
  xi <- xiFunRatioMean(m=m, nx, ny, x, y)
  
  nu <- which(xi > u)
  if (length(nu) > 0) {
    mStop <- m[min(nu)]
  } else {
    mStop <- mMax
  }

  if (plot) {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(x = m, y = xi, pch=19, ylab = expression(paste(xi, "(m)")), 
         main= bquote(paste(m["stop"]^"asym", " = ", .(mStop), sep = "")),
         cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4)
    abline(a = u, b = 0, col="red")
  }
  
  class(mStop) <- "mStop"
  return(mStop)
}

xiFunDiffMean <- function(m, nx, ny, x, y){
  #' Utitlity function for mStop
  #'
  #' Calculates xi(m)
  #' @export
  
  N <- nx + ny
  # lambda <- nx / N
  
  xbar <- mean(x)
  ybar <- mean(y)
  
  Emean <- gLin(m = m, nx, ny, xbar = xbar, ybar = ybar)
    
  Evar <-  VarRLin(m, nx, ny, x, y)
  
  t <- xbar - ybar

  xi <- (t - Emean) / sqrt(Evar)
  
  return(xi)
}

mStopDiffMean <- function(x, y, B=1000, plot = FALSE){
  #' Estimate partition at which fastPerm will stop for difference statistics
  #'
  #' This function estimates the partition at which fastPerm will
  #' will stop, based on asymptotic approximation.
  #' @param x First sample
  #' @param y Second sample
  #' @param B Number of Monte Carlo iterations within each partition. Defaults to 1,000.
  #' @keywords fast permutation stopping partition
  #' @export
  #' @examples
  #' x <- rnorm(100, 1)
  #' y <- rexp(100, 0)
  #' mStopDiffMean(x, y)

  # swap x and y if xbar < ybar
  if (mean(x) < mean(y)) {
    ytemp <- y
    y <- x
    x <- ytemp
  }

  xbar <- mean(x)
  ybar <- mean(y)
  
  nx <- length(x)
  ny <- length(y)
  
  pmf <- PiLgamma(nx, ny)
  mMax <- min(which(pmf == max(pmf)))

  u <- qnorm(1-1/B)
  m <- 1:min(mMax,(min(nx, ny)-1))
  
  xi <- xiFunDiffMean(m=m, nx, ny, x, y)
  
  nu <- which(xi > u)
  if (length(nu) > 0) {
    mStop <- m[min(nu)]
  } else {
    mStop <- mMax
  }
  
  if (plot) {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(x = m, y = xi, pch=19, ylab = expression(paste(xi, "(m)")), 
         main= bquote(paste(m["stop"]^"asym", " = ", .(mStop), sep = "")),
         cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4)
    abline(a = u, b = 0, col="red")
  }
  
  class(mStop) <- "mStop"
  return(mStop)
}

print.mStop <- function(mStopObj){
  #' Print function for mStop functions
  #'
  #' This function prints the expected stopping partition
  #' @param mStopObj Output from the mStop functions
  #' @keywords mStop print
  #' @export
  #' @examples
  #' x <- rnorm(110, 1)
  #' y <- rnorm(95, 0)
  #' print(mStopDiffMean(x, y))
  
  result <- paste("Expected mStop = ", mStopObj, sep = "")
  writeLines(result)
}
