#' An algorithm for least-squares curve fitting.
#'
#' This algorithm provides a numerical solution to the problem of maximizing a
#' function. This is more efficient than the Gauss-Newton-like algorithm when
#' starting from points very far from the final maximum. A new convergence test
#' is implemented (RDM) in addition to the usual stopping criterion : stopping
#' rule is when the gradients are small enough in the parameters metric
#' (GH-1G).
#'
#' Convergence criteria are very strict as they are based on derivatives of the
#' log-likelihood in addition to the parameter and log-likelihood stability. In
#' some cases, the program may not converge and reach the maximum number of
#' iterations fixed at 500.  In this case, the user should check that parameter
#' estimates at the last iteration are not on the boundaries of the parameter
#' space.  If the parameters are on the boundaries of the parameter space, the
#' identifiability of the model should be assessed.  If not, the program should
#' be run again with other initial values, with a higher maximum number of
#' iterations or less strict convergence tolerances.
#'
#' @param b an optional vector containing the initial values for the
#' parameters. Default is 0.1 for every parameter.
#' @param m an optional parameter if the vector of parameter is not missing
#' compulsory if b is not given.
#' @param fn The function to be maximized, with argument the
#' vector of parameters over which maximization is to take place.  It should
#' return a scalar result.
#' @param maxiter optional maximum number of iterations for the marqLevAlg
#' iterative algorithm. Default is 500.
#' @param epsa optional threshold for the convergence criterion based on the
#' parameter stability. Default is 0.001.
#' @param epsb optional threshold for the convergence criterion based on the
#' log-likelihood stability. Default is 0.001.
#' @param epsd optional threshold for the relative distance to maximum. This
#' criterion has the nice interpretation of estimating the ratio of the
#' approximation error over the statistical error, thus it can be used for
#' stopping the iterative process whathever the problem. Default is 0.01.
#' @param digits Number of digits to print in outputs. Default value is 8.
#' @param print.info Logical.Equals to TRUE if report (parameters at iteration,
#' function value, convergence criterion ...) at each iteration is requested.
#' Default value is FALSE.
#' @param blinding Logical. Equals to TRUE if the algorithm is allowed to go on
#' in case of an infinite or not definite value of function. Default value is
#' FALSE.
#' @param multipleTry Integer, different from 1 if the algorithm is allowed to
#' go for the first iteration in case of an infinite or not definite value of
#' gradients or hessian. This account for a starting point to far from the
#' definition set. As many tries as requested in multipleTry will be done by
#' changing the starting point of the algorithm. Default value is 25.
#' @param nproc number of processors for parallel computing
#' @param clustertype one of the supported types from \code{\link[parallel]{makeCluster}}
#' @param .packages for parallel setting only, packages used in the fn function
#' @param \dots other arguments of the fn function 
#'
#' @return \item{cl}{ summary of the call to the function marqLevAlg.  }
#' \item{ni}{ number of marqLevAlg iterations before reaching stopping
#' criterion.  } \item{istop}{ status of convergence: =1 if the convergence
#' criteria were satisfied, =2 if the maximum number of iterations was reached,
#' =4 if the algorithm encountered a problem in the function computation.  }
#' \item{v}{ vector containing the upper triangle matrix of variance-covariance
#' estimates at the stopping point.  } \item{fn.value}{ function evaluation at
#' the stopping point.  } \item{b}{ stopping point value.  } \item{ca}{
#' convergence criteria for parameters stabilisation.  } \item{cb}{ convergence
#' criteria for function stabilisation.  } \item{rdm}{ convergence criteria on
#' the relative distance to minimum.  } \item{time}{ a running time.  }
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
#'
#' @export
#'
#' @examples
#'
#'
#' ### 1
#' ### initial values
#' b <- c(8,9)
#' ### your function
#' f1 <- function(b){
#' 	return(-4*(b[1]-5)^2-(b[2]-6)^2)
#' }
#' ## Call
#' test1 <- marqLevAlg(b=b,fn=f1)
#'
#'microbenchmark::microbenchmark(marqLevAlg(b=b,fn=f1),
#'                               marqLevAlg(b=b,fn=f1,nproc = 3)
#'                               )
#'
#'
#'
#' ### 2
#' ## initial values
#' b <- c(3,-1,0,1)
#' ## your function
#' f2 <- function(b){
#' 	return(-((b[1]+10*b[2])^2+5*(b[3]-b[4])^2+(b[2]-2*b[3])^4+10*(b[1]-b[4])^4))
#' }
#' ## Call
#' test2 <- marqLevAlg(b=b,fn=f2)
#' test2$b
#'
#' test2_par <- marqLevAlg(b=b,fn=f2, nproc = 2)
#' test2_par$b
#'microbenchmark::microbenchmark(marqLevAlg(b=b,fn=f2),
#'                               marqLevAlg(b=b,fn=f2, nproc = 3)
#'                               )
#'


marqLevAlg <- function(b,m=FALSE,fn,gr=NULL,hess=NULL,maxiter=500,epsa=0.001,epsb=0.001,epsd=0.01,digits=8,print.info=FALSE,blinding=TRUE,multipleTry=25,nproc=1,clustertype=NULL,file="",.packages,...){
	cl <- match.call()
	if (missing(m) & missing(b)) stop("The 'marqLevAlg' alogorithm needs a vector of parameters 'b' or his length 'm'")
	if(missing(m)) m <- length(b)	
	if(missing(b)) b <- rep(0.1,m)
	if(length(b) != m){
		if(length(b) < m){
			b.temp <-NULL
			b.temp <- c(b.temp,b,rep(b[1],(m-length(b))))
			b <- b.temp 
		}else{
			m <- length(b)
		}
	}  

	if(missing(fn)) stop("The argument 'funcpa' is missing.")

        
        if(nproc>1)
            {
                if(is.null(clustertype))
                {
                    clustpar <- parallel::makeCluster(nproc)
                }
                else
                {
                    clustpar <- parallel::makeCluster(nproc, type=clustertype)
                }
                
                doParallel::registerDoParallel(clustpar)
            }


	funcpa <- function(b,...){fn(b,...)} #modif -fn en fn
	if(!is.null(gr)) grad <- function(b,...){gr(b,...)}
	if(!is.null(hess)) hessian <- function(b,...){hess(b,...)}

	flush.console()
	ptm <- proc.time()
	#cat("\n")
	#cat("Be patient. The program is computing ...\n")
	
###initialisation
	binit <- b
	th <- 1e-5
	eps <- 1e-7
	nfmax <- m*(m+1)/2    
	ca <- epsa+1
	cb <- epsb+1
	rl1 <- -1.e+10    
	ni <- 0
	istop <- 0
	da <- 1E-2
	dm <- as.double(5)
	nql <- 1
	m1 <- m*(m+1)/2
	ep <- 1E-20
	delta <- rep(0,m)
	b1 <- rep(0,m)
	idpos <- 0
	v <- rep(0,m*(m+3)/2)
	ind.func1 <- 0
	fu <- rep(0,(m*(m+3)/2))
	gonflencountmax <- 10

## old parameters iteration -1
	old.b <- b
	old.rl <- 0
	old.ca <- 1
	old.cb <- 1
	old.dd <- 1	
## 	
	repeat{	
	
		if (sum(!is.finite(b))>0){

			cat("Probably too much accuracy requested...\n")
			cat("Last step values :\n")
			cat("      b :",round(old.b,digits),"\n")
			cat("      likelihood :",round(-old.rl,digits),"\n")
			cat("      Convergence criteria: parameters stability=", round(old.ca,digits), "\n")
			cat("                          : likelihood stability=", round(old.cb,digits), "\n") 
			cat("                          : best relative distance to maximum obtained (RDM)=", round(old.dd,digits), "\n")
			stop("")
			 
		}
		res.out.error <- list("old.b"=round(old.b,digits),"old.rl"=round(old.rl,digits),"old.ca"=round(old.ca,digits),"old.cb"=round(old.cb,digits),"old.dd"=round(old.dd,digits))
	
		if(is.null(gr)){ 
			deriv <- deriva(nproc,b,funcpa,.packages=.packages,...)
			v <- deriv$v
			rl <- deriv$rl
			
			if((multipleTry > 1) & (ni ==0)){
				kk <- 0
				while(((kk < multipleTry) & (!is.finite(rl)))){
					kk <- kk + 1
					b <- b/2
					deriv <- deriva(nproc,b,funcpa,.packages=.packages,...)
					v <- deriv$v
					rl <- deriv$rl
				}
			} 
		}else{
			v <- NULL
			rl=funcpa(b,...)
			
			if(is.null(hess)){
				deriv <- deriva_grad(b,grad,...)
				v <- c(v,deriv$hessian,grad(b,...))
			}else{
				tmp.hessian <- hessian(b) 
				if(is.matrix(tmp.hessian)){
					tmp.hessian <- tmp.hessian[upper.tri(tmp.hessian,diag=T)]
				}
				v <- c(v,tmp.hessian,grad(b,...))	
			}
			
		}
		if((sum(is.finite(b))==m) && !is.finite(rl)){
			cat("Problem of computation. Verify your function specification...\n")
			cat("Infinite likelihood with finite parameters : b=",round(old.b,digits),"\n")
			cat("      - Check the computation and the continuity,\n")
			cat("      - Check that you minimize the function.\n")
			stop("")
		
		}

		if(((sum(!is.finite(b)) > 0) || (sum(!is.finite(rl)) > 0)) && (ni==0)){
			cat("Problem of computation. Verify your function specification...\n")
			cat("First check b length in parameters specification.\n")
			stop("")
		}
		rl1 <- rl      
		dd <- 0 
			
		for(i in 1:m){
			for(j in i:m){
				ij <- (j-1)*j/2+i
				fu[ij]=v[ij]
			}
		}
		fu.int <- fu[1:(m*(m+1)/2)]
	
		dsinv <- .Fortran("dsinv",fu.out=as.double(fu.int),as.integer(m),as.double(ep),ier=as.integer(0),det=as.double(0),PACKAGE="marqLevAlgParallel")
				
		ier <- dsinv$ier
		fu[1:(m*(m+1)/2)] <- dsinv$fu.out
		if (ier == -1){
			dd <- epsd+1 
			v_tmp <- v[(m*(m+1)/2+1):(m*(m+3)/2)]
			#dd <- sum(v_tmp*v_tmp)
		}else{
			dd <- ghg(m,v,fu)$ghg/m
                        if(is.na(dd)) dd <- epsd+1 ##??? changmt
		}
		
        if(print.info){
		cat("------------------ iteration ",ni,"------------------\n",file=file,append=TRUE)
		cat("Log_likelihood ",round(rl,digits),"\n",file=file,append=TRUE)
		cat("Convergence criteria: parameters stability=", round(ca,digits), "\n",file=file,append=TRUE)
		cat("                    : likelihood stability=", round(cb,digits), "\n",file=file,append=TRUE) 
		## if (ier == -1){
		## 	cat("                    : Matrix inversion for RDM failed \n",file=file,append=TRUE)	
		## }else{
		## 	cat("                    : Matrix inversion for RDM successful \n",file=file,append=TRUE)
		## }
		cat("                    : relative distance to maximum(RDM)=", round(dd,digits), "\n",file=file,append=TRUE)

		nom.par <- paste("parameter",c(1:m),sep="")
		id <- 1:m
		indice <- rep(id*(id+1)/2)
		Var <- fu[indice]
		SE <- sqrt(abs(Var))
		res.info <- data.frame("coef"=round(b,digits),"SE.coef"=round(SE,digits),"Var.coef"=round(Var,digits))
		rownames(res.info) <- nom.par
                if(file=="") print(res.info)
                else write.table(res.info,file=file,append=TRUE)
		cat("\n")
	}
		old.b <- b
		old.rl <- rl
		old.ca <- ca
		old.cb <- cb
                #if(is.na(dd)) browser()
		if(dd<=old.dd){old.dd <- dd}
		if((ca < epsa) & (cb < epsb) & (dd < epsd)){break}

		
		tr <- 0
		for(i in 1:m){
			ii <- i*(i+1)/2
			tr <- tr+abs(v[ii])
		}
		tr <- tr/m
		
		ncount <- 0
		ga <- 0.01

		fu <- v
		
		for(i in 1:m){
			ii <- i*(i+1)/2
			if (v[ii] != 0){
				fu[ii] <- v[ii]+da*((1.e0-ga)*abs(v[ii])+ga*tr)
			}else{
				fu[ii] <- da*ga*tr
			}
		}
			
		dchole <- .Fortran("dchole",fu=as.double(fu),as.integer(m),as.integer(nql),idpos=as.integer(0),PACKAGE="marqLevAlgParallel")
		fu <- dchole$fu
		idpos <- dchole$idpos

		while(idpos != 0){
		
			ncount <- ncount + 1
			if((ncount <= 3) | (ga >= 1)){
				da <- da * dm
			}else{
				ga <- ga * dm
				if(ga > 1) ga <- 1
			}
			
			fu <- v
			
			for(i in 1:m){
				ii <- i*(i+1)/2
				if (v[ii] != 0){
					fu[ii] <- v[ii]+da*((1.e0-ga)*abs(v[ii])+ga*tr)
				}else{
					fu[ii] <- da*ga*tr
				}
			}
			
		dchole <- .Fortran("dchole",fu=as.double(fu),as.integer(m),as.integer(nql),idpos=as.integer(0),PACKAGE="marqLevAlgParallel")
		
			idpos <- dchole$idpos
			fu <- dchole$fu
			if (ncount >= gonflencountmax){
				break
			} 

		}
		delta <- fu[(nfmax+1):(nfmax+m)]
		b1 <- b + delta
		rl <- funcpa(b1,...)
		
		if(blinding){
			if(is.na(rl)){
				cat("rl :",rl,"\n")
				rl <- -500000
			}
		}else{
			if(is.na(rl)){
				cat(" Probably wrong definition of function FN \n")
				cat("      ->  invalid number (infinite or NA)\n")
				cat("          value of function is :",round(-rl,digits),"\n")
				
				istop <- 4
				break
			}
		}
		if (rl1 < rl){
		
			if(da < eps){
				da <- eps
			}else{
				da <- da/(dm+2)
			}
			goto800 <- func1(b,rl1,rl,delta,ni,maxiter)
			ca <- goto800$ca
			cb <- goto800$cb
			b <- goto800$b
			ni <- goto800$ni
			ind.func1 <- 1
			if(ni >= maxiter){
				istop <- 2
				break
			}
		}else{
			maxt <- max(abs(delta)) 
	
			if(maxt == 0){
				vw <- th
			}else{
				vw <- th/maxt
			}
		
			step <- log(1.5)

			sears <- searpas(vw,step,b,delta,funcpa,res.out.error,...)

			fi <- sears$fi
			vw <- sears$vw
			rl <- -fi
			if(rl == -1.e9){
				istop <- 4
				break
			}
			delta <- vw*delta
			da <- (dm-3)*da
			goto800 <- func1(b,rl1,rl,delta,ni,maxiter)
 			ca <- goto800$ca
			cb <- goto800$cb
			b <- goto800$b
			ni <- goto800$ni
			   
			if(ni >= maxiter){
				istop <- 2
				break
			}
		}

	}
        
        if(nproc>1)
        {
            parallel::stopCluster(clustpar)
        }

        
	if((istop %in% 2:4)==F) istop <- 1
	cost <- proc.time() - ptm
	result <- list(cl=cl,ni=ni,ier=ier,istop=istop,v=fu[1:(m*(m+1)/2)],fn.value=rl,b=b,ca=ca,cb=cb,rdm=dd,time=round(cost[3],3)) 
	class(result) <- "marqLevAlg"
        if(print.info==TRUE)
            {
                if(file!="")
                    {
                        fileres <- sub(".txt","_last.txt",file)
                        dput(result,file=fileres)
                    }
            }
	

	return(result)
}























