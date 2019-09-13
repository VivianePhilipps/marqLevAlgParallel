\name{deriva}
\alias{deriva}

\title{
Numerical approch to derivate.
}
\description{
The function to return the first, second derivate and the information score matrix. There are the central finite-difference and forward finite-difference will be used.

}
\usage{
deriva(nproc=1, b, funcpa, .packages=NULL, ...)
}

\arguments{
\item{nproc}{number of processors for parallel computing
}
  \item{b}{value of parameters to be optimized over.
}
  \item{funcpa}{
function to be minimized (or maximized), with argument the vector of parameters over which minimization is
to take place.  It should return a scalar result.
}
\item{.packages}{character vector of packages that funcpa depends on}
\item{\dots}{other arguments of the funcpa function}
}

\value{
  \item{v}{
    the information score matrix.
  }
  \item{rl}{
    log-likelihood or likelihood of the model.
  }
}
\references{
Donald W. Marquardt An algorithm for Least-Squares Estimation of Nonlinear Parameters. Journal of the Society for Industrial and Applied Mathematics, Vol. 11, No. 2. (Jun, 1963), pp. 431-441.
}
\author{
Viviane Philipps, Boris Hejblum, Cecile Proust-Lima, Daniel Commenges
}

\examples{

b <- 0.1
f <- function(b){return((2*b[1]**2+3*b[1]))}

d <- deriva(b=b,funcpa=f)

}


