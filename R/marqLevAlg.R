#' A parallelized general-purpose optimization based on Marquardt-Levenberg algorithm
#'
#' This algorithm provides a numerical solution to the problem of unconstrained local
#' optimization. This is more efficient than the Gauss-Newton-like algorithm when
#' starting from points very far from the final maximum. A new convergence test
#' is implemented (RDM) in addition to the usual stopping criterion : stopping
#' rule is when the gradients are small enough in the parameters metric
#' (GH-1G).
#'
#' Convergence criteria are very strict as they are based on derivatives of the
#' objective function in addition to the parameter and objective function stability. In
#' some cases, the program may not converge and reach the maximum number of
#' iterations fixed at 500.  In this case, the user should check that parameter
#' estimates at the last iteration are not on the boundaries of the parameter
#' space.  If the parameters are on the boundaries of the parameter space, the
#' identifiability of the model should be assessed.  If not, the program should
#' be run again with other initial values, with a higher maximum number of
#' iterations or less strict convergence tolerances. An alternative is to remove some
#' parameters from the Hessian matrix.
#'
#' @param b an optional vector containing the initial values for the
#' parameters. Default is 0.1 for every parameter.
#' @param m number of parameters. Optional if b is specified.
#' @param fn the function to be optimized, with first argument the
#' vector of parameters over which optimization is to take place (argument b).
#' It should return a scalar result.
#' @param gr a function to return the gradient value for a specific point.
#' If missing, finite-difference approximation will be used.
#' @param hess a function to return the hessian matrix for a specific point.
#' If missing, finite-difference approximation will be used.
#' @param maxiter optional maximum number of iterations for the marqLevAlg
#' iterative algorithm. Default is 500.
#' @param epsa optional threshold for the convergence criterion based on the
#' parameter stability. Default is 0.0001.
#' @param epsb optional threshold for the convergence criterion based on the
#' objective function stability. Default is 0.0001.
#' @param epsd optional threshold for the relative distance to maximum. This
#' criterion has the nice interpretation of estimating the ratio of the
#' approximation error over the statistical error, thus it can be used for
#' stopping the iterative process whathever the problem. Default is 0.0001.
#' @param partialH optional vector giving the indexes of the parameters to be dropped from
#' the Hessian matrix to define the relative distance to maximum/minimum. If specified,
#' this option will only be considered at iterations where the two first convergence
#' criteria are satisfied (epsa and epsb) and if the total Hessian is not invertible.
#' By default, no partial Hessian is defined.
#' @param digits number of digits to print in outputs. Default value is 8.
#' @param print.info logical indicating if information about computation should be
#' reported at each iteration.
#' Default value is FALSE.
#' @param blinding logical. Equals to TRUE if the algorithm is allowed to go on
#' in case of an infinite or not definite value of function. Default value is
#' FALSE.
#' @param multipleTry integer, different from 1 if the algorithm is allowed to
#' go for the first iteration in case of an infinite or not definite value of
#' gradients or hessian. As many tries as requested in multipleTry will be done by
#' changing the starting point of the algorithm. Default value is 25.
#' @param nproc number of processors for parallel computing
#' @param clustertype one of the supported types from \code{\link[parallel]{makeCluster}}
#' @param file optional character giving the name of the file where the outputs
#' of each iteration should be written (if print.info=TRUE).
#' @param .packages for parallel setting only, packages used in the fn function
#' @param minimize logical indicating if the fn function should be minimized or maximized. By default minimize=TRUE, the function is minimized.
#' @param \dots other arguments of the fn, gr and hess functions 
#'
#' @return \item{cl}{ summary of the call to the function marqLevAlg.  }
#' \item{ni}{ number of marqLevAlg iterations before reaching stopping
#' criterion.  } \item{istop}{ status of convergence: =1 if the convergence
#' criteria were satisfied, =2 if the maximum number of iterations was reached,
#' =3 if convergence criteria with partial Hessian matrix were satisfied,
#' =4 if the algorithm encountered a problem in the function computation.  }
#' \item{v}{if istop=1 or istop=3, vector containing the upper triangle matrix of variance-covariance
#' estimates at the stopping point. Otherwise v contains the second derivatives
#' of the fn function with respect to the parameters.} \item{grad}{vector containing the gradient
#' at the stopping point.} \item{fn.value}{ function evaluation at
#' the stopping point.  } \item{b}{ stopping point value.  } \item{ca}{
#' convergence criteria for parameters stabilisation.  } \item{cb}{ convergence
#' criteria for function stabilisation.  } \item{rdm}{ convergence criteria on
#' the relative distance to minimum (or maximum).  } \item{time}{ a running time.  }
#' @author Melanie Prague, Viviane Philipps, Cecile Proust-Lima, Boris Hejblum, Daniel Commenges, Amadou Diakite
#' @references \emph{marqLevAlg Algorithm}
#'
#' Donald W. marquardt An algorithm for Least-Squares Estimation of Nonlinear
#' Parameters. Journal of the Society for Industrial and Applied Mathematics,
#' Vol. 11, No. 2. (Jun, 1963), pp. 431-441.
#'
#' \emph{Convergence criteria : Relative distance to Minimim (or Maximum)}
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
#' ### example 1
#' ### initial values
#' b <- c(8,9)
#' ### your function
#' f1 <- function(b){
#' 	return(-4*(b[1]-5)^2-(b[2]-6)^2)
#' }
#' ### gradient
#' g1 <- function(b){
#'      return(c(-8*(b[1]-5),-2*(b[2]-6)))
#' }
#' ## Call
#' test1 <- mla(b=b, fn=f1, minimize=FALSE)
#'
#' \dontrun{
#'microbenchmark::microbenchmark(mla(b=b, fn=f1, minimize=FALSE),
#'                               mla(b=b, fn=f1, minimize=FALSE, nproc=2),
#'                               mla(b=b, fn=f1, gr=g1, minimize=FALSE),
#'                               mla(b=b, fn=f1, gr=g1, minimize=FALSE, nproc=2),
#'                               times=10)
#'         }
#'
#'
#'
#' ### example 2
#' ## initial values
#' b <- c(3,-1,0,1)
#' ## your function
#' f2 <- function(b){
#' 	return(-((b[1]+10*b[2])^2+5*(b[3]-b[4])^2+(b[2]-2*b[3])^4+10*(b[1]-b[4])^4))
#' }
#' ## Call
#' test2 <- mla(b=b, fn=f2, minimize=FALSE)
#' test2$b
#'
#' test2_par <- mla(b=b, fn=f2, minimize=FALSE, nproc=2)
#' test2_par$b
#' 
#' \dontrun{
#'microbenchmark::microbenchmark(mla(b=b, fn=f2, minimize=FALSE),
#'                               mla(b=b, fn=f2, minimize=FALSE, nproc=2),
#'                               times=10)
#'         }
#'
#'
#'\dontrun{
#'### example 3 : a linear mixed model
#'## the log-likelihood is implemented in the loglikLMM function
#'## the gradient is implemented in the gradLMM function
#'
#'## data
#'Y <- dataEx$Y
#'X <- as.matrix(cbind(1,dataEx[,c("t","X1","X3")],dataEx$t*dataEx$X1))
#'ni <- as.numeric(table(dataEx$i))
#'
#'## initial values
#'binit <- c(0,0,0,0,0,1,1)
#'
#'## estimation in sequential mode, with numeric derivatives
#'estim <- marqLevAlg(b=binit, fn=loglikLMM, minimize=FALSE, X=X, Y=Y, ni=ni)
#'## estimation in parallel mode, with numeric derivatives
#'estim2 <- marqLevAlg(b=binit, fn=loglikLMM, minimize=FALSE, X=X, Y=Y, ni=ni, 
#'nproc=2, clustertype="FORK")
#'## estimation in sequential mode, with analytic gradient
#'estim3 <- marqLevAlg(b=binit, fn=loglikLMM, gr=gradLMM, minimize=FALSE, X=X, Y=Y, ni=ni)
#'## estimation in parallel mode, with analytic gradient
#'estim4 <- marqLevAlg(b=binit, fn=loglikLMM, gr=gradLMM, minimize=FALSE, X=X, Y=Y, ni=ni, 
#'nproc=2, clustertype="FORK")
#'}

marqLevAlg <- function(b,m=FALSE,fn,gr=NULL,hess=NULL,maxiter=500,epsa=0.0001,epsb=0.0001,epsd=0.0001,partialH=NULL,digits=8,print.info=FALSE,blinding=TRUE,multipleTry=25,nproc=1,clustertype=NULL,file="",.packages=NULL,minimize=TRUE,...){
	cl <- match.call()
	if (missing(m) & missing(b)) stop("The 'marqLevAlg' algorithm needs a vector of parameters 'b' or his length 'm'")
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

        
        if(nproc>1){
            if(is.null(clustertype)){
                clustpar <- parallel::makeCluster(nproc)#, outfile="")
            }
            else{
                clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
            }
            
            doParallel::registerDoParallel(clustpar)
        }


        if(minimize==TRUE){
             funcpa <- function(b,...){-fn(b,...)} 
            if(!is.null(gr)) grad <- function(b,...){-gr(b,...)}            
        }
        else{
            funcpa <- function(b,...){fn(b,...)} 
            if(!is.null(gr)) grad <- function(b,...){gr(b,...)}
        }
	if(!is.null(hess)) hessian <- function(b,...){hess(b,...)}

	flush.console()
	ptm <- proc.time()
	
###initialisation
	binit <- b
	th <- 1e-5
	eps <- 1e-7
	nfmax <- m*(m+1)/2    
	ca <- epsa+1
	cb <- epsb+1
        dd <- epsd+1
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
        ier <- 0
## 	
	repeat{	
	
		if (sum(!is.finite(b))>0){

                    cat("Infinite parameters...\n")
                    cat("Last step values :\n")
                    cat("      b :",round(old.b,digits),"\n")
                    if(minimize){cat("      objective function :",round(-old.rl,digits),"\n")
                    } else {cat("      objective function :",round(old.rl,digits),"\n")}
                    cat("      Convergence criteria: parameters stability=", round(old.ca,digits), "\n")
                    cat("                          : function stability=", round(old.cb,digits), "\n") 
                    cat("                          : best relative distance to maximum obtained (RDM)=", round(old.dd,digits), "\n")
                    istop <- 4
                    rl <- -1.e9
                    break
			 
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
				deriv <- deriva_grad(nproc=nproc,b,grad,.packages=.packages,...)
				v <- c(v,deriv$hessian,grad(b,...))
			}else{
				tmp.hessian <- hessian(b,...) 
				if(is.matrix(tmp.hessian)){
					tmp.hessian <- tmp.hessian[upper.tri(tmp.hessian,diag=T)]
				}
				v <- c(v,tmp.hessian,grad(b,...))	
			}
			
		}
		if((sum(is.finite(b))==m) && !is.finite(rl)){
			cat("Problem of computation. Verify your function specification...\n")
			cat("Infinite value with finite parameters : b=",round(old.b,digits),"\n")
                        istop <- 4
                        rl <- -1.e9
                        break
		}

		if(((sum(!is.finite(b)) > 0) || (sum(!is.finite(rl)) > 0)) && (ni==0)){
			cat("Problem of computation. Verify your function specification...\n")
			cat("Infinite value or parameters\n")
			istop <- 4
                        rl <- -1.e9
                        break
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
	
		dsinv <- .Fortran(C_dsinv,fu.out=as.double(fu.int),as.integer(m),as.double(ep),ier=as.integer(0),det=as.double(0))
				
                ier <- dsinv$ier
                fu[1:(m*(m+1)/2)] <- dsinv$fu.out
		if (ier == -1){
			dd <- epsd+1 
		}else{
			dd <- ghg(m,v,fu)$ghg/m
                        if(is.na(dd)) dd <- epsd+1
		}
		
        if(print.info){
		cat("------------------ iteration ",ni,"------------------\n",file=file,append=TRUE)
		cat("Function value ",round(rl,digits),"\n",file=file,append=TRUE)
		cat("Convergence criteria: parameters stability=", round(ca,digits), "\n",file=file,append=TRUE)
		cat("                    : function stability=", round(cb,digits), "\n",file=file,append=TRUE) 
		cat("                    : relative distance to maximum(RDM)=", round(dd,digits), "\n",file=file,append=TRUE)

		nom.par <- paste("parameter",c(1:m),sep="")
		res.info <- data.frame("coef"=round(b,digits))
		rownames(res.info) <- nom.par
                if(file=="") print(res.info)
                else write.table(res.info,file=file,append=TRUE)
		cat("\n")
	}
		old.b <- b
		old.rl <- rl
		old.ca <- ca
		old.cb <- cb
		if(dd<=old.dd){old.dd <- dd}
		if((ca < epsa) & (cb < epsb) & (dd < epsd)){break}

            if((ca < epsa) & (cb < epsb) & (ier == -1) & (sum(partialH)>0)){

                ## partial Hessian matrix 
                mr <- m - length(partialH)
                H <- matrix(NA, m, m)
                H[upper.tri(H, diag=TRUE)] <- v[1:(m*(m+1)/2)]
                Hr <- H[setdiff(1:m, partialH),setdiff(1:m, partialH)]
                fur <- Hr[upper.tri(Hr, diag=TRUE)]
                
                ## inversion of the partial Hessian matrix
                dsinvr <- .Fortran(C_dsinv,fu.out=as.double(fur),as.integer(mr),as.double(ep),ier=as.integer(0),det=as.double(0))

                
                ier <- dsinvr$ier
		fur[1:(mr*(mr+1)/2)] <- dsinvr$fu.out
		if (ier == -1){
                    dd <- epsd+1
                    if(print.info){
                        cat("Inversion of partial Hessian matrix failed \n \n")
                    }
		}else{
                    ## third convergence criteria with partial H
                    vr <- rep(NA, mr*(mr+3)/2)
                    ir <- 0
                    for(i in 1:m) {
                        if(!(i %in% partialH)){
                            ir <- ir+1
                            vr[mr*(mr+1)/2+ir] <- v[m*(m+1)/2+i]
                        }
                    }
                    dd <- ghg(mr,vr,fur)$ghg/mr
                    if(is.na(dd)) dd <- epsd+1
		}

                if(print.info){
                    cat("RDM with partial Hessian matrix= ", dd, "\n \n")
                }
                
                if(dd < epsd) {
                    ## convergence with partial H
                    istop <- 3
                    v <- rep(NA, m*(m+3)/2)
                    v[m*(m+1)/2 + setdiff(1:m, partialH)] <- vr[mr*(mr+1)/2+1:mr]
                    Hr <- matrix(NA, mr, mr)
                    Hr[upper.tri(Hr, diag=TRUE)] <- fur
                    H <- matrix(NA, m, m)
                    H[setdiff(1:m, partialH), setdiff(1:m, partialH)] <- Hr
                    fu <- H[upper.tri(H, diag=TRUE)]
                    break
                }
                
            }
            
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
			
		dchole <- .Fortran(C_dchole,fu=as.double(fu),as.integer(m),as.integer(nql),idpos=as.integer(0))
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
			
		dchole <- .Fortran(C_dchole,fu=as.double(fu),as.integer(m),as.integer(nql),idpos=as.integer(0))
		
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
				if(minimize)  cat("rl :",-rl,"\n") else cat("rl :",rl,"\n")
				rl <- -500000
			}
		}else{
			if(is.na(rl)){
				cat(" Numerical problem by computing fn \n")
				if(minimize) {cat("          value of function is :",round(-rl,digits),"\n")
                                }else{
                                    cat("          value of function is :",round(rl,digits),"\n")
				}
				istop <- 4
                                rl <- -1.e9
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
        
        if(nproc>1){
            parallel::stopCluster(clustpar)
        }
        
        if(minimize==TRUE) rl <- -rl
	if((istop %in% 2:4)==F) istop <- 1
	cost <- proc.time() - ptm

        result <- list(cl=cl,ni=ni,ier=ier,istop=istop,v=fu[1:(m*(m+1)/2)],
                       grad=v[(m*(m+1)/2)+1:m],fn.value=rl,b=b,
                       ca=ca,cb=cb,rdm=dd,time=round(cost[3],3))
	class(result) <- "marqLevAlg"
        
        if(print.info==TRUE){
            if(file!=""){
                fileres <- sub(".txt","_last.txt",file)
                dput(result,file=fileres)
            }
        }
	

	return(result)
}





#' @rdname marqLevAlg
#' @export
mla <- marqLevAlg


















