ipf2 <-
function(rt=NULL,ct=NULL,m=matrix(1,length(rt),length(ct)),tol=1e-05,maxit=500,iter=FALSE){
  if(!is.null(rt) & !is.null(ct))
    if(round(sum(rt))!=round(sum(ct))) 
      stop("row and column totals are not equal, ensure sum(rt)==sum(ct)")
  n<-list(i=rt,
          j=ct)
  mu<-m
  mu.marg<-n
  m.fact<-n
  if(iter==TRUE)
    rd<-paste("%.",nchar(format(tol,scientific=FALSE))-2,"f",sep="")
  it<-0; max.diff<-tol*2
  while(it==0 | max.diff>tol & it<maxit ){
    if(!is.null(ct)){
      mu.marg$j <- apply(mu,2,sum)
      m.fact$j <- n$j/mu.marg$j
      m.fact$j[is.nan(m.fact$j)]<-0
      mu <- sweep(mu, 2, m.fact$j, "*")
    }
    if(!is.null(rt)){
      mu.marg$i <- apply(mu,1,sum)
      m.fact$i <- n$i/mu.marg$i
      m.fact$i[is.nan(m.fact$i)]<-0
      mu <- sweep(mu, 1, m.fact$i, "*")
    }
    it<-it+1
    max.diff<-max(abs(unlist(n)-unlist(mu.marg)))
    if(iter==TRUE)
      cat(sprintf(rd,unlist(m.fact)), fill = T)
  }
  return(list(mu=mu,it=it,tol=max(abs(unlist(n)-unlist(mu.marg)))))
}
