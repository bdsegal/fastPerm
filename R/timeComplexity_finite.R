# Functions derived in the Appendix -----------------------

muw <- function(m, lambda, N, xbar, ybar) {
	#' Utitlity function for mStop
	#'
	#' Calculates expected value of W in partition m
	
  m / sqrt(lambda * (1 - lambda) * N) * (xbar-ybar)
}

g <- function(m, nx, ny, xbar, ybar){
	#' Utitlity function for mStop
	#'
	#' Calculates g function of the expected value of W
	
	N <- nx + ny
	lambda <- nx / N
	
	numerator <- sqrt((lambda * N)/(1 - lambda) ) * xbar + 
                     muw(m, lambda, N, xbar, ybar)
	
	denominator <- sqrt((1 - lambda) * N / lambda) * ybar - 
                      muw(m, lambda, N, xbar, ybar)
	
	return((1 - lambda) / lambda * numerator / denominator)
}

gp <- function(m, nx, ny, xbar, ybar){
	#' Utitlity function for mStop
	#'
	#' Calculates derivative of the g function of the expected value of W
	
	N <- nx + ny
	lambda <- nx / N
	
	numerator <- sqrt((lambda * N) / (1 - lambda)) * xbar + 
               sqrt(((1 - lambda) * N) / lambda) * ybar

	denominator <- (sqrt(((1 - lambda) * N) / lambda) * ybar - 
                  muw(m, lambda, N, xbar, ybar))^2
	
	return((1 - lambda) / lambda * numerator / denominator)
}

ER <- function(m, nx, ny, xbar, ybar){
	#' Utitlity function for mStop
	#'
	#' Calculates expected value of the ratio stat in partition m
	
	num <- xbar + m / nx * (ybar - xbar)
	den <- ybar + m / ny * (xbar - ybar)
	
	return(num / den)
}
	
VarW <- function(m, nx, ny, x, y){
	#' Utitlity function for mStop
	#'
	#' Calculates variance of W
	
	N <- nx + ny
	lambda <- nx / N
	
	# variance terms
	EAy <- 1 / (lambda * (1 - lambda) * N) * m / ny * (1 - m / ny) *
	        sum(y^2)
	EAx <- 1 / (lambda * (1 - lambda) * N) * m / nx * (1 - m / nx) *
          sum(x^2)
					
	# get cross product terms
	yCrossProd <- outer(y,y)
	yCrossProdSum <- sum(yCrossProd) - sum(diag(yCrossProd))
	
	xCrossProd <- outer(x,x)
	xCrossProdSum <- sum(xCrossProd) - sum(diag(xCrossProd))
	
	# covariance terms
	ECy <- 1 / (lambda * (1 - lambda) * N) * (m * (ny - m)) / 
          (ny^2 * (ny - 1)) * yCrossProdSum
			
	ECx <- 1 / (lambda * (1 - lambda) * N) * (m * (nx - m)) /
          (nx^2 * (nx - 1)) * xCrossProdSum

	return(EAx + EAy - ECx - ECy)
}

VarR <- function(m, nx, ny, x, y, xbar, ybar){
	#' Utitlity function for mStop
	#'
	#' Calculates variance of the ratio stat R
	
	gp(m, nx, ny, xbar, ybar)^2 * 
  VarW(m, nx, ny, x, y)
}

xiFunRatioMean <- function(m, nx, ny, x, y){
	#' Utitlity function for mStop
	#'
	#' Calculates xi(m)
	
	N <- nx + ny
	lambda <- nx / N
	
	xbar <- mean(x)
	ybar <- mean(y)
	
	Emean <- g(m, nx, ny, xbar, ybar)
	Evar <- VarR(m, nx, ny, x, y, xbar, ybar)
	
	t <- xbar / ybar

	xi <- (t - Emean) / sqrt(Evar)
	
	return(xi)
}

mStopRatioMean <- function(x, y, B=1000, plot = TRUE){
	#' Estimate partition at which fastPerm will stop for ratio statistics
	#'
	#' This function estimates the partition at which fastPerm will
	#' will stop, based on asymptotic approximation.
	#' @param x First sample
	#' @param y Second sample
	#' @param B Number of Monte Carlo iterations within each partition. Defaults to 1000.
	#' @keywords fast permtution stopping partition
	#' @export
	#' @examples
	#' x <- rexp(100, 4)
	#' y <- rexp(100, 2)
	#' mStopRatioMean(x, y)

	# swap x and y if ybar < xbar
	if (mean(y) < mean(x)) {
		ytemp <- y
		y <- x
		x <- ytemp
	}

	nx <- length(x)
	ny <- length(y)

	u <- qnorm(1-1/B)
	m <- 1:(nx-1)
	xi <- xiFunRatioMean(m=m, nx, ny, x, y)
	
	mStop <- m[min(which(xi > u))]
	
	if (plot) {
		plot(x = m, y = xi, pch=19, ylab = expression(paste(xi, "(m)")), 
		main= paste("Expected mStop = ", mStop, sep = ""))
		abline(a = u, b = 0, col="red")
	}
	
	class(mStop) <- "mStop"
	return(mStop)
}

xiFunDiffMean <- function(m, nx, ny, x, y){
	#' Utitlity function for mStop
	#'
	#' Calculates xi(m)
	
	N <- nx + ny
	lambda <- nx / N
	
	xbar <- mean(x)
	ybar <- mean(y)
	
	Emean <- xbar - ybar + 1 / sqrt(lambda * (1 - lambda) * N) *
		muw(m, lambda, N, xbar, ybar)
		
	Evar <- 1 / (lambda * (1 - lambda) * N) * VarW(m, nx, ny, x, y)
	
	t <- xbar - ybar

	xi <- (t - Emean) / sqrt(Evar)
	
	return(xi)
}

mStopDiffMean <- function(x, y, B=1000, plot = TRUE){
	#' Estimate partition at which fastPerm will stop for difference statistics
	#'
	#' This function estimates the partition at which fastPerm will
	#' will stop, based on asymptotic approximation.
	#' @param x First sample
	#' @param y Second sample
	#' @param B Number of Monte Carlo iterations within each partition. Defaults to 1000.
	#' @keywords fast permtution stopping partition
	#' @export
	#' @examples
	#' x <- rnorm(100, 1)
	#' y <- rexp(100, 0)
	#' mStopDiffMean(x, y)

	# swap x and y if ybar < xbar
	if (mean(y) < mean(x)) {
		ytemp <- y
		y <- x
		x <- ytemp

	}

	xbar <- mean(x)
	ybar <- mean(y)
	
	nx <- length(x)
	ny <- length(y)
	
	u <- qnorm(1-1/B)
	m <- 1:(nx-1)
	
	xi <- xiFunDiffMean(m=m, nx, ny, x, y)
	
	mStop <- m[min(which(xi > u))]
	
	if (plot) {
		plot(x = m, y = xi, pch=19, ylab = expression(paste(xi, "(m)")), 
		main= paste("Expected mStop = ", mStop, sep = ""))
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
	#' x <- rnorm(100, 1)
	#' y <- rnorm(100, 0)
	#' print(mStopDiffMean(x, y))
	
		result <- paste("Expected mStop = ", mStopObj, sep = "")
		writeLines(result)
}