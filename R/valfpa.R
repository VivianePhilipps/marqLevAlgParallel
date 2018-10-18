valfpa <- function(vw,b,delta,funcpa,...){
    if(is.na(vw)) return(-2E9)
	bk <- b + (exp(vw)*delta)
        #cat("dans valfpa, b= ",b,"\n")
        #cat("  vw= ",vw,"   delta=",delta,"\n")
        #cat("             bk= ",bk,"\n")
	return(-funcpa(bk,...)) 
}
