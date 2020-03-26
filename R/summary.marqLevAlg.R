#' @export
summary.marqLevAlg <- function(object,digits=8,...){
    x <- object
    if (!inherits(x, "marqLevAlg")) stop("use only with \"marqLevAlg\" objects")
    
    
    cat(" \n")
    cat("                   Robust marqLevAlg algorithm                   ", "\n")
    cat(" \n")
    cl <- x$cl
    dput(cl)
    cat(" \n")
    cat("Iteration process:", "\n")
    cat("      Number of parameters:", length(x$b)," \n")
    cat("      Number of iterations:", x$ni, "\n")
    cat("      Optimized objective function:", round(x$fn.value,digits)," \n")
    if(x$istop==1) cat("      Convergence criteria satisfied","\n")
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
    cat("                    : relative distance to maximum(RDM)=", round(x$rdm,digits), "\n")
    if(x$istop!=4&x$istop!=5) {
        cat(" \n")
        cat("Final parameter values:", "\n")
        id <- 1:length(x$b)
        indice <- rep(id*(id+1)/2)
        se <-sqrt(x$v[indice])
        wald <- (x$b/se)**2
        z <- abs(qnorm((1 + .95)/2))
        binf <- x$b-1.96*se
        bsup <- x$b+1.96*se
        
        tmp <- data.frame("coef"=format(round(x$b,3)),"SE coef"=format(round(se,3)),"Wald"=format(wald,4),"P-value"=round(1 - pchisq(wald, 1),5),"binf"=round(binf,3),"bsup"=round(bsup,3))
        print(tmp,row.names=F)
        cat(" \n")
    }
}