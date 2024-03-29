% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deriva_grad.R
\name{deriva_grad}
\alias{deriva_grad}
\title{Numerical derivatives of the gradient function}
\usage{
deriva_grad(nproc = 1, b, grad, .packages = NULL, .export = NULL, ...)
}
\arguments{
\item{nproc}{number of processors for parallel computing}

\item{b}{value of parameters to be optimized over}

\item{grad}{the gradient of the function to be minimized (or maximized)}

\item{.packages}{character vector of packages that grad depends on}

\item{.export}{character vector of objects/functions that grad depends on}

\item{\dots}{other arguments of the grad function}
}
\value{
\item{hessian}{vector containing the upper part of the information score matrix}
}
\description{
The function computes the information score matrix in the case where the first
derivatives of the function to optimize are analytically known. Therefore,
minus the derivatives of the gradient are computed by central finite differences.
}
\author{
Viviane Philipps, Boris Hejblum, Cecile Proust-Lima, Daniel Commenges
}
