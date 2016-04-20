

fastPermAsym <- function(x, y, testStat = ratioMean, plot = FALSE){
#' Fast approximation of small permutation p-values via asymptotic results
#'
#' This function approximates the p-value (two-sided) for a two
#' sample permutation test by partitioning the permutation space,
#' and using asymptotic results on the p-values within each partition.
#' This function is most useful for small p-values, e.g. p < 10^-6, and
#' large sample sizes.
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
#' fastPerm(x, y)
#'
#' x <- rnorm(110, 0, 1)
#' y <- rnorm(110, 1, 1)
#' mStopDiffMean(x, y)
#' fastPerm(x, y, testStat = diffMean)
  
  comparison <- attributes(testStat)$comparison  
  if (comparison == "ratio" & (min(x,y)<0)) {
    stop("Some sample values <=0. With ratio statistics, all values must be greater than 0.")
  }
  
  if(comparison == "ratio"){
    xiFun <- xiFunRatioMean
  } else if (comparison == "difference"){
    xiFun <- xiFunDiffMean
  } else {
    stop("testStat must be either ratioMean or diffMean")
  }

  # swap x and y if xbar < ybar
  if (mean(x) < mean(y)) {
    ytemp <- y
    y <- x
    x <- ytemp
  }
  
  nx <- length(x)
  ny <- length(y)
  nMax <- min(nx, ny)
  
  t0 <- testStat(x,y)
  
  pmf <- PiLgamma(nx, ny)
  mMax <- min(which(pmf == max(pmf)))

  if (length(mMax)==1 & mMax[1]==nMax){
    m <- 1:mMax
  } else if (length(mMax)==1 & mMax[1]<nMax){
    m <- c(1:mMax,((mMax-1):0)[1:(nMax-mMax)])
  } else if (length(mMax)==2 & mMax[2]==nMax){
    m <- c(1:mMax[1], mMax[1])
  } else {
    m <- c(1:mMax[1], (mMax[1]:0)[1:(nMax - mMax[1])])
  }
  
  xi <- xiFun(m=m, nx, ny, x, y)
 
  # if nx != ny, set p-value in 0 partition to 1 and 
  # estimate p-value in partition min(nx,ny)
  if (nx != ny) { 
    pNorm <- c(1, pnorm(xi, lower.tail=FALSE)) %*% pmf
    pT <- c(1, pt(xi, df=nx+ny, lower.tail=FALSE)) %*% pmf
    
  # if nx==ny, set both the 0 and nx partition to 1
  } else { 
    pNorm <- c(1, pnorm(xi[-length(xi)], lower.tail=FALSE), 1) %*% pmf
    pT <- c(1, pt(xi[-length(xi)], df=nx+ny, lower.tail=FALSE), 1) %*% pmf
  }

  
#   if (plot) {
#     plot(x = m, y = xi, pch=19, ylab = expression(paste(xi, "(m)")), 
#     main= paste("Expected mStop = ", mStop, sep = ""))
#     abline(a = u, b = 0, col="red")
#   }
  fpAsym <- list(pNorm=as.numeric(pNorm),
                 pT=as.numeric(pT),
                 comparison = attributes(testStat)$comparison,
                 summary = attributes(testStat)$summary,
                 t0=t0)

  class(fpAsym) <- "fastPermAsym"
  return(fpAsym)
}

print.fastPermAsym <- function(fp){
  #' Print function for fastPerm
  #'
  #' This function prints the results of fastPermAsym
  #' @param fp Output from the fastPermAsym function
  #' @keywords fastPerm print
  #' @export
  #' @examples
  #' x <- rexp(100, 5)
  #' y <- rexp(100, 2)
  #' print(fastPermAsym(x, y, testStat = ratioMean))
  
  result <- paste("    fastPerm Two Sample Test   \n\n",
                  fp$comparison, " of ", fp$summary, "s\n", 
                  "\nobserved statistic = ", signif(fp$t0,3),
                  "\n\np-value \t Approximation", 
                  "\n-------------------------------",
                  "\n", signif(fp$pT,3), "\t t-distribution",
                  "\n", signif(fp$pNorm,3), "\t normal",
                  sep = "")
  
  writeLines(result)
}
