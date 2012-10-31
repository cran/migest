cm2 <-
function(rt=NULL,ct=NULL,m=matrix(1,length(rt),length(ct)),tol=1e-05,maxit=500,iter=TRUE)
{
  if(round(sum(rt))!=round(sum(ct))) 
    stop("row and column totals are not equal, ensure sum(rt)==sum(ct)")
  i<-dim(m)[1];  j<-dim(m)[2]
  alpha <- rep(1,i)
  beta <- rep(1,j)
  if(iter==TRUE){
    rd<-paste("%.",nchar(format(tol,scientific=FALSE))-2,"f",sep="")
    cat(sprintf(rd,c(alpha,beta)), fill = T)
  }
  alpha.old <- alpha+1; beta.old <- beta+1
  it<-1;  max.diff<-tol*2
  while(max.diff>tol & it<maxit ){
    beta.old <- beta
    for(j in 1:j) {
      beta[j] <- ct[j]/sum(alpha * m[, j])
    }
    alpha.old <- alpha
    for(i in 1:i) {
      alpha[i] <- rt[i]/sum(beta * m[i,  ])
    }
    it<-it+1
    max.diff<-max(abs(alpha-alpha.old), abs(beta-beta.old))
    if(iter==TRUE)
      cat(sprintf(rd,c(alpha,beta)), fill = T)
  }
  return(list(N=alpha%*%t(beta)*m,
              theta=c(mu=1,alpha=alpha,beta=beta)))
}
