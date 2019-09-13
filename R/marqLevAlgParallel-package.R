#' A parallelized algorithm for least-squares curve fitting.
#'
#' This algorithm provides a numerical solution to the problem of maximizing a
#' function. This is more efficient than the Gauss-Newton-like algorithm when
#' starting from points very far from the final maximum. A new convergence test
#' is implemented (RDM) in addition to the usual stopping criterion : stopping
#' rule is when the gradients are small enough in the parameters metric
#' (GH-1G).
#'
#' \tabular{ll}{ Package: \tab marqLevAlgParallel\cr Type: \tab Package\cr
#' Version: \tab 1.0\cr Date: \tab 2018-10-18\cr License: \tab GPL (>= 2.0)\cr
#' LazyLoad: \tab yes\cr } This algorithm provides a numerical solution to the
#' problem of maximizing a function. This is more efficient than the
#' Gauss-Newton-like algorithm when starting from points very far from the final
#' maximum. A new convergence test is implemented (RDM) in addition to the
#' usual stopping criterion : stopping rule is when the gradients are small
#' enough in the parameters metric (GH-1G).
#'
#' @name marqLevAlgParallel-package
#' @docType package
#' @author Viviane Philipps, Cecile Proust-Lima, Boris Hejblum
#' @references \emph{marqLevAlg Algorithm}
#'
#' Donald W. marquardt An algorithm for Least-Squares Estimation of Nonlinear
#' Parameters. Journal of the Society for Industrial and Applied Mathematics,
#' Vol. 11, No. 2. (Jun, 1963), pp. 431-441.
#'
#' \emph{Convergence criteria : Relative distance to Minimim}
#'
#' Commenges D. Jacqmin-Gadda H. Proust C. Guedj J. A Newton-like algorithm for
#' likelihood maximization the robust-variance scoring algorithm
#' arxiv:math/0610402v2 (2006)
#' @keywords marqLevAlg algorithm optimization maximisation
#' package
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom stats pchisq
#' @importFrom stats qnorm
#' @importFrom utils flush.console
#' @importFrom utils write.table
#' @useDynLib marqLevAlgParallel
NULL
