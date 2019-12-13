\name{weib}
\alias{weib}
\docType{data}
\title{
Simulated dataset for the weibull function.
}
\description{
a dataset contains 936 rows and 5 columns.
}

\usage{data(weib)}
\format{
  A data frame with 936 observations on the following 5 variables.
  \describe{
    \item{\code{c}}{a numeric vector.}
    \item{\code{t0}}{entry time.}
    \item{\code{t1}}{left border of interval censored.}
    \item{\code{t2}}{right border of interval censored.}
    \item{\code{ve}}{a numeric vector.}
  }
}

\examples{
data(weib)
head(weib)
}
\keyword{datasets}
