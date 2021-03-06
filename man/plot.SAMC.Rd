% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAMCfunctions.R
\name{plot.SAMC}
\alias{plot.SAMC}
\title{plot function for SAMC}
\usage{
\method{plot}{SAMC}(sam)
}
\arguments{
\item{sam}{Output from the SAMC function}
}
\description{
This function plots the convergence diagnostics for the SAMC algorithm and prints
out the associated values. If the SAMC converged, then it should have sampled
nearly uniformly from all regions. The maximum discrepancy is the maximum difference
between the number of times SAMC sampled from a region and the expected amount if the
sampling were uniform. Yu et al. suggest that the maximum discrepancy should be less 
than 0.2.
}
\examples{
x <- rexp(100, 5)
y <- rexp(100, 2)
sam <- SAMC(x, y)
plot(sam)
}
\keyword{SAMC}
\keyword{convergence}
\keyword{diagnostics}
\keyword{summary}

