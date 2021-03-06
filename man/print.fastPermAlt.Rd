% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastPermAltFunctions.R
\name{print.fastPermAlt}
\alias{print.fastPermAlt}
\title{Print function for fastPerm}
\usage{
\method{print}{fastPermAlt}(fp)
}
\arguments{
\item{fp}{Output from the fastPerm function}
}
\description{
This function prints the results of fastPerm
}
\examples{
x <- rexp(100, 5)
y <- rexp(100, 2)
print(fastPerm(x, y, testStat = ratioMean))
}
\keyword{fastPerm}
\keyword{print}

