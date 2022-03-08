#' Summary of optimization
#'
#' A short summary of parameters estimates by marqLevAlg algorithm.
#'
#' @param object a marqLevAlg object
#' @param digits Number of digits to print in outputs. Default value is 8.
#' @param loglik Logical indicating if the objective function is a log-likelihood. By
#' default, loglik=FALSE.
#' @param \dots other (unsued) arguments
#'
#' @return A data frame containing as many rows as estimated parameters. If
#'  loglik=FALSE, it includes one column containing the estimated
#'  parameters values. If loglik=TRUE, it includes 6 columns : the
#'  estimated parameters, their standard errors, the corresponding Wald
#'  statistic, the associated p-value and the boundaries of the 95\% confidence
#'  interval.
#'
#' @seealso \code{\link{marqLevAlg}}, \code{\link{print.marqLevAlg}}
#'
#' @keywords summary
#'
#' @author V. Philipps, C. Proust-Lima, B. Hejblum, D. Commenges, M. Prague, A. Diakite
#'
#' @examples
#'
#' f1 <- function(b){	
#'	return(4*(b[1]-5)^2+(b[2]-6)^2)	
#' }
#' test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,
#' epsd=0.001,fn=f1)
#'
#' summary(test.marq)
#' 
#' @export
summary.marqLevAlg <- function(object,digits=8,loglik=FALSE,...){
    x <- object
    if (!inherits(x, "marqLevAlg")) stop("use only with \"marqLevAlg\" objects")
    
    
    cat(" \n")
    cat("                   Robust marqLevAlg algorithm                   ", "\n")
    cat(" \n")
    cl <- x$cl
    minimize <- TRUE
    if(length(cl$minimize)){
        if(cl$minimize==FALSE) minimize <- FALSE
    }
    dput(cl)
    cat(" \n")
    cat("Iteration process:", "\n")
    cat("      Number of parameters:", length(x$b)," \n")
    cat("      Number of iterations:", x$ni, "\n")
    cat("      Optimized objective function:", round(x$fn.value,digits)," \n")
    if(x$istop==1) cat("      Convergence criteria satisfied","\n")
    if(x$istop==3) cat("      Convergence criteria with partial Hessian matrix satisfied","\n")
    if(x$istop==2) cat("      Maximum number of iteration reached without convergence","\n")
    if(x$istop==4|x$istop==5)  {
        cat("      The program stopped abnormally. No results can be displayed.\n")
    }
    cat(" \n")
    cat("Convergence criteria: parameters stability=", round(x$ca[1],digits), "\n")
    cat("                    : objective function stability=", round(x$cb,digits), "\n") 
    if (x$ier == -1){
        cat("                    : Matrix inversion for RDM failed \n")	
    }else{
        cat("                    : Matrix inversion for RDM successful \n")
    }
    if(minimize==TRUE){
        cat("                    : relative distance to minimum(RDM)=", round(x$rdm,digits), "\n")
    }else{
        cat("                    : relative distance to maximum(RDM)=", round(x$rdm,digits), "\n")
    }

    if(x$istop!=4&x$istop!=5) {
        cat(" \n")
        cat("Final parameter values:", "\n")
        id <- 1:length(x$b)
        indice <- rep(id*(id+1)/2)
        se <- sqrt(x$v[indice])
        wald <- (x$b/se)**2
        z <- abs(qnorm((1 + .95)/2))
        binf <- x$b-1.96*se
        bsup <- x$b+1.96*se

        if(loglik==FALSE){
            tmp <- data.frame("coef"=format(round(x$b,3)))
            cat(format(round(x$b,3)),"\n")
        }else{
            tmp <- data.frame("coef"=format(round(x$b,3)),"SE coef"=format(round(se,3)),"Wald"=format(wald,4),"P-value"=round(1 - pchisq(wald, 1),5),"binf"=round(binf,3),"bsup"=round(bsup,3))
            print(tmp,row.names=F)
        }
        cat(" \n")

        return(invisible(tmp))
    }
}
