#' @export
summary.marqLevAlg <- function(object,digits=8,...)
{
    x <- object
    if (!inherits(x, "marqLevAlg")) stop("use only with \"marqLevAlg\" objects")

    cat(" \n")
    cat("Values:", "\n")
    cat(" \n")
    id <- 1:length(x$b)
    indice <- rep(id*(id+1)/2)
    se <-sqrt(x$v[indice])
    wald <- (x$b/se)**2
    z <- abs(qnorm((1 + .95)/2))
    binf <- x$b-1.96*se
    bsup <- x$b+1.96*se
    tmp <- data.frame("coef"=format(round(x$b,3)),"SE coef"=format(round(se,3)),"Wald"=format(wald,4),"P-value"=round(1 - pchisq(wald, 1),5),"binf"=round(binf,3),"bsup"=round(bsup,3))
   # print(tmp,row.names=FALSE)
    
    cat(" \n")
    cat("Number of iterations: ", x$ni, "\n")
    cat(" \n")
    cat("Convergence criteria: parameters stability=", round(x$ca[1],digits), "\n")
    cat("                    : likelihood stability=", round(x$cb,digits), "\n") 
    cat("                    : relative distance to maximum=", round(x$rdm,digits), "\n")
    cat(" \n")
    cat("Goodness-of-fit statistics:", "\n")
    cat("      maximum log-likelihood:", round(x$fn.value,digits)," \n")
    cat(" \n")
    cat(" \n")
    cat(" \n")
}

