\name{summary.marqLevAlg}
\alias{summary.marqLevAlg}

\title{
summary of optimization.
}
\description{
A short summary of parameters estimates by marqLevAlg algorithm.
}

\usage{
\method{summary}{marqLevAlg}(object, digits, loglik=FALSE, \dots)
}

\arguments{
  \item{object}{
a marqLevAlg object.
}
  \item{digits}{
Number of digits to print in outputs. Default value is 8.
}
  \item{loglik}{
Logical indicating if the objective function is a log-likelihood. By
  default, loglik=FALSE.
}
  \item{\dots}{
other unusued arguments.
}
}

\value{A data frame containing as many rows as estimated parameters. If
  loglik=FALSE, it includes one column containing the estimated
  parameters values. If loglik=TRUE, it includes 6 columns : the
  estimated parameters, their standard errors, the corresponding Wald
  statistic, the associated p-value and the boundaries of the 95\% confidence
  interval.}

\seealso{
 \code{\link{marqLevAlg}}, 
 \code{\link{print.marqLevAlg}}
}
  
\author{
V. Philipps, C. Proust-Lima, B. Hejblum, D. Commenges, M. Prague, A. Diakite
}

\examples{
f1 <- function(b){	
	return(4*(b[1]-5)^2+(b[2]-6)^2)	
}
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,
epsd=0.001,fn=f1)

summary(test.marq)
}

\keyword{summary}
