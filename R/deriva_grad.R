

deriva_grad <- function(nproc=1,b,grad,.packages=NULL,...){
    m <- length(b)
    bh2 <- b
    bh <- b 
    v <- rep(0,m*(m+1)/2)
    k <- 0
    if(nproc>1)
    {
        ## remplacer les 2 boucles par une seule
        grid <- cbind(rep(1:m,1:m),unlist(sapply(1:m,function(k) seq(1,k))))
        mm <- nrow(grid)
        h <- sapply(b,function(x){max(1E-7,(1E-4*abs(x)))})

        ## derivees du gradient
        v <- foreach(k=1:(m*(m+1)/2),
                      .combine=c,
                      .export=c("grid","b","h","grad"),
                      .packages=.packages) %dopar%
            {
                i <- grid[k,2]
                j <- grid[k,1]
                
                bp <- b
                bp[j] <- b[j]+h[j]
                av <- grad(bp,...)[i]
                
                bm <- b
                bm[j] <- b[j]-h[j]
                ar <- grad(bm,...)[i]
                
                d <- (ar-av)/(2*h[j])               
                d
            }
    }
    else
    {
	for(j in 1:m){
            for(i in 1:j){
                k <- k+1
                bh <- b 
                bh2 <- b
                thj <- -max(1e-7, 1e-4 * abs(b[j]))
                bh2[j] <- b[j]-thj
                bh[j] <- b[j]+thj
                v[k] <- (grad(bh,...)[i]-grad(bh2,...)[i])/(2*(-thj))
            }
	}
    }
    
    result <- list(hessian=v)
    return(result)
}

