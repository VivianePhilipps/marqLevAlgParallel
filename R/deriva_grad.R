

deriva_grad <- function(nproc=1,b,grad,.packages=NULL,...){
    m <- length(b)
    if(nproc>1)
    {
        h <- sapply(b,function(x){max(1E-7,(1E-4*abs(x)))})

        ## derivees du gradient
        vtmp <- foreach(k=1:m,
                        .combine=rbind,
                        .export=c("grid","b","h","grad"),
                        .packages=.packages) %dopar%
            {
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
        h <- sapply(b,function(x){max(1E-7,(1E-4*abs(x)))})
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

