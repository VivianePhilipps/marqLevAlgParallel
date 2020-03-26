#' Log-likelihood of a linear mixed model with random intercept
#' 
#' @param b numeric vector specifying the parameter's values in the
#' following order : first the fixed effects and then the standard
#' deviation of the random intercept and of the independent error
#' @param Y numeric vector including the dependent outcome vector ordered by subject
#' @param X numeric matrix including the covariates
#' @param ni interger vector giving the number of repeated measures for each subject
#' @return the log-likelihood of the linear mixed model at point b
#' @export

loglikLMM <- function(b,Y,X,ni)
{
  beta <- b[1:ncol(X)]
  alpha <- b[ncol(X)+1]
  sigma <- b[ncol(X)+2]
  
  loglik <- 0
  j <- 0
  
  for(i in 1:length(ni))
  {
    Vi <- alpha^2 * matrix(1,nrow=ni[i],ncol=ni[i]) + sigma^2 * diag(ni[i])
    
    loglik <- loglik - (ni[i]*log(2*pi) + log(det(Vi)) + 
        t(Y[j+1:ni[i]] - X[j+1:ni[i],]%*%beta) %*% solve(Vi) %*% (Y[j+1:ni[i]] - X[j+1:ni[i],]%*%beta) )/2 
    
    j <- j + ni[i]
  }
  
  return(as.numeric(loglik))
}


#' Gradient of the log-likelihood of a linear mixed model with random intercept
#' 
#' @param b numeric vector specifying the parameter's values in the
#' following order : first the fixed effects and then the standard
#' deviation of the random intercept and of the independent error
#' @param Y numeric vector including the dependent outcome vector ordered by subject
#' @param X numeric matrix including the covariates
#' @param ni interger vector giving the number of repeated measures for each subject
#' @return a vector containing the gradient of the log-likelihood of the linear mixed model at point b
#' @export
gradLMM <- function(b,Y,X,ni)
{
  beta <- b[1:ncol(X)]
  alpha <- b[ncol(X)+1]
  sigma <- b[ncol(X)+2]
  
  gradbeta <- 0
  gradalpha <- 0
  gradsigma <- 0
  j <- 0
  
  for(i in 1:length(ni))
  {
    Vi <- alpha^2 * matrix(1,nrow=ni[i],ncol=ni[i]) + sigma^2 * diag(ni[i])
    invVi <- solve(Vi)
    
    gradbeta <- gradbeta + t(X[j+1:ni[i],,drop=FALSE]) %*% invVi %*% (Y[j+1:ni[i]] - X[j+1:ni[i],,drop=FALSE]%*%beta)
    
    gradalpha <- gradalpha -( sum(diag(2*alpha*invVi%*%matrix(1,nrow=ni[i],ncol=ni[i]))) + 
        t(Y[j+1:ni[i]] - X[j+1:ni[i],,drop=FALSE]%*%beta) %*% (-2*alpha* invVi %*% matrix(1,nrow=ni[i],ncol=ni[i]) %*% invVi) %*% (Y[j+1:ni[i]] - X[j+1:ni[i],,drop=FALSE]%*%beta))/2

    gradsigma <- gradsigma -( sum(diag(2*sigma*invVi)) + 
        t(Y[j+1:ni[i]] - X[j+1:ni[i],,drop=FALSE]%*%beta) %*% (-2*sigma* invVi %*% invVi) %*% (Y[j+1:ni[i]] - X[j+1:ni[i],,drop=FALSE]%*%beta))/2

    j <- j + ni[i]
  }
  
  return(c(gradbeta,gradalpha,gradsigma))
}


