ipf3 <-
function(rt=NULL,ct=NULL,m=NULL,tol=1e-05,maxit=500,iter=TRUE){
  if(any(round(colSums(rt))!=round(rowSums(ct))))
    stop("row and column totals are not equal for one or more sub-tables, ensure colSums(rt)==rowSums(ct)")
  n<-list(ik=rt,
          jk=t(ct))
  R<-dim(rt)[1]
  
  #set up offset
  if(is.null(m)){
    m<-array(1,c(dim(rt),dim(rt)[1]))
    dimnames(m)<-list(orig=dimnames(rt)[[1]],dest=dimnames(ct)[[1]],pob=dimnames(rt)[[2]])
  }

  mu<-m
  mu.marg<-n
  m.fact<-n
  it<-0; max.diff<-tol*2
  while(max.diff>tol & it<maxit){
    mu.marg$ik <- apply(mu,c(1,3),sum)
    m.fact$ik <- n$ik/mu.marg$ik
    m.fact$ik[is.nan(m.fact$ik)]<-0
    m.fact$ik[m.fact$ik==Inf]<-0
    mu <- sweep(mu, c(1,3), m.fact$ik, "*")
    
    mu.marg$jk <- apply(mu, c(2,3), sum)
    m.fact$jk <- n$jk/mu.marg$jk
    m.fact$jk[is.nan(m.fact$jk)]<-0
    mu <- sweep(mu, c(2,3), m.fact$jk, "*")
    
    it<-it+1
    max.diff<-max(abs(unlist(n)-unlist(mu.marg)))
    if(iter==TRUE)
      cat(c(it, max.diff), "\n")
  }
  return(list(mu=mu,it=it,tol=max.diff))
}
