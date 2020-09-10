#' A parallelized general-purpose optimization based on Marquardt-Levenberg algorithm
#'
#' This algorithm provides a numerical solution to the problem of
#' miniminzing/maximizing a function. This is more efficient than the
#' Gauss-Newton-like algorithm when starting from points very far from the final
#' minimum/maximum. A new convergence test
#' is implemented (RDM) in addition to the usual stopping criterion : stopping
#' rule is when the gradients are small enough in the parameters metric
#' (GH^{-1}G).
#'
#' 
#' \Sexpr[stage=build,results=hide]{descr <- packageDescription("marqLevAlg")}
#' 
#' \tabular{ll}{ Package: \tab marqLevAlg\cr Type: \tab Package\cr
#' Version: \tab \Sexpr[stage=build]{descr$Version} \cr
#' Date: \tab \Sexpr[stage=build]{descr$Date} \cr License: \tab GPL (>= 2.0)\cr
#' LazyLoad: \tab yes\cr } This algorithm provides a numerical solution to the
#' problem of optimizing a function. This is more efficient than the
#' Gauss-Newton-like algorithm when starting from points very far from the final
#' maximum. A new convergence test is implemented (RDM) in addition to the
#' usual stopping criterion : stopping rule is when the gradients are small
#' enough in the parameters metric (GH-1G).
#'
#' @name marqLevAlg-package
#' @docType package
#' @author Viviane Philipps, Cecile Proust-Lima, Boris Hejblum, Melanie Prague, Daniel Commenges, Amadou Diakite
#' @references \emph{marqLevAlg Algorithm}
#'
#' Philipps V. Hejblum B.P. Prague M. Commenge D. Proust-Lima C.
#' Robust and Efficient Optimization Using a Marquardt-Levenberg Algorithm
#' with R Package marqLevAlg
#'
#' Donald W. marquardt An algorithm for Least-Squares Estimation of Nonlinear
#' Parameters. Journal of the Society for Industrial and Applied Mathematics,
#' Vol. 11, No. 2. (Jun, 1963), pp. 431-441.
#'
#' \emph{Convergence criteria : Relative distance to Maximum}
#'
#' Commenges D. Jacqmin-Gadda H. Proust C. Guedj J. A Newton-like algorithm for
#' likelihood maximization : the robust-variance scoring algorithm
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
#' @useDynLib marqLevAlg, .registration=TRUE, .fixes="C_"
NULL
