#' Numerical derivatives of the gradient function
#'
#' The function computes the information score matrix in the case where the first
#' derivatives of the function to optimize are analytically known. Therefore,
#' minus the derivatives of the gradient are computed by central finite differences.
#'
#' @param nproc number of processors for parallel computing
#' @param b value of parameters to be optimized over
#' @param grad the gradient of the function to be minimized (or maximized)
#' @param .packages character vector of packages that grad depends on
#' @param \dots other arguments of the grad function
#'
#' @return \item{hessian}{vector containing the upper part of the information score matrix}
#' @author Viviane Philipps, Boris Hejblum, Cecile Proust-Lima, Daniel Commenges
#'
#' @export
#' 
deriva_grad <- function(nproc=1,b,grad,.packages=NULL,...){
    m <- length(b)
    h <- sapply(b,function(x){max(1E-7,(1E-4*abs(x)))})
    if(nproc>1)
    {
        ## derivees du gradient
        vtmp <- foreach(j=1:m,
                        .combine=rbind,
                        .packages=.packages) %dopar%
            {
                k <- j
                bp <- b
                bp[k] <- bp[k] + h[k]
                av <- grad(bp,...)
                
                bm <- b
                bm[k] <- bm[k] - h[k]
                ar <- grad(bm,...)
                
                d <- (ar-av)/(2*h[k])               
                d
            }
    }
    else
    {
        vtmp <- matrix(NA,m,m)
	for(j in 1:m){
            bp <- b
            bp[j] <- bp[j] + h[j]
            bm <- b
            bm[j] <- bm[j] - h[j]
            vtmp[j,] <- (grad(bm,...)-grad(bp,...))/(2*h[j])
	}
    }

    v <- vtmp[upper.tri(vtmp,diag=TRUE)]
    result <- list(hessian=v)
    return(result)
}

