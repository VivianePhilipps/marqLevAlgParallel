% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.marqLevAlg.R
\name{print.marqLevAlg}
\alias{print.marqLevAlg}
\title{Summary of a \code{marqLevAlg} object}
\usage{
\method{print}{marqLevAlg}(x, digits = 8, ...)
}
\arguments{
\item{x}{a marqLevAlg object.}

\item{digits}{Number of digits to print in outputs. Default value is 8.}

\item{\dots}{other (unusued) arguments.}
}
\description{
The function provides a summary of a \code{marqLevAlg} optimisation.
}
\examples{

f1 <- function(b){	
return(4*(b[1]-5)^2+(b[2]-6)^2)	
}
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,
epsd=0.001,fn=f1)

test.marq

}
\seealso{
\code{link{summary.marqLevAlg}}
}
\author{
V. Philipps, C. Proust-Lima, B. Hejblum, D. Commenges, M. Prague, A. Diakite
}
\keyword{print}
