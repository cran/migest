ffs <-
function(P1,P2,d,b,m=NULL,net=FALSE,d.adj=FALSE,...){
  R<-nrow(P1)
  #set up offset
  if(is.null(m)){
    m<-array(1,c(dim(P1),dim(P1)[1]))
  }
  if(is.null(dimnames(m))){
    dimnames(m)<-list(orig=dimnames(P1)[[1]],dest=dimnames(P2)[[1]],pob=dimnames(P1)[[2]])
  }
  if(net==TRUE){
    if(round(sum(P2-P1),1)!=round(sum(b-d),1)){
      message("sum(P2-P1): ", sum(P2-P1))
      message("sum(b-d):   ", sum(b-d))
      stop("difference in stock tables must be equal to sum of all births - sum of all deaths")
    }
  } 
  
  #set up migration matrix with births/deaths and others rows and columns
  y<-m
  y<-array(0,dim(m)+c(2,2,0))
  dimnames(y)<-list(o=c(dimnames(m)[[1]],"b","O"),d=c(dimnames(m)[[1]],"d","O"),b=dimnames(m)[[3]])
  y[,,]<-0
  
  #step 1-2a take off deaths
  d.mat<-ipf2(ctot=d,m=P1)$mu
  if(d.adj==TRUE){
    if(net==TRUE){
      stop("set net to FALSE, currently have two conflicting methods to balance data")
    }
    d.mat<-ipf2(ctot=d,rtot=-1*(rowSums(P2)-(rowSums(P1)+b)),m=P1)$mu  
  }
  y[1:R,R+1,]<-t(d.mat)
  P1.adj<-P1-d.mat
  
  #step 1-2b take off births
  b.mat<-diag(b)
  y[R+1,1:R,]<-b.mat
  P2.adj<-P2-b.mat
  
  #step 3-4a take off moves in from external or adjust P1.adj rows
  dif<-rowSums(P1.adj) - rowSums(P2.adj)
  if(net==FALSE){
    in.mat<-t(ipf2(ctot=pmax(dif,0),m=t(P2.adj))$mu)
    P1.adj<-P1.adj-in.mat
    y[R+2,1:R,]<-t(in.mat)
  }
  if(net==TRUE){
    P1.adj<-ipf2(rtot=rowSums(P1.adj)-dif/2,ctot=colSums(P1.adj),m=P1.adj)$mu
  }
  
  #step 3-4b take off moves out from external or adjust P2.adj rows
  if(net==FALSE){
    out.mat<-t(ipf2(ctot=pmax(-dif,0),m=t(P1.adj))$mu)
    P2.adj<-P2.adj-out.mat
    y[1:R,R+2,]<-t(out.mat)
  }
  if(net==TRUE){
    P2.adj<-ipf2(rtot=rowSums(P2.adj)+dif/2,ctot=colSums(P2.adj),m=P2.adj)$mu
  }
  
  #step 5 calculate
  #ipf<-ipf3.qi(rtot=t(P1.adj),ctot=P2.adj,m=m)
  ipf<-ipf3.qi(rtot=t(P1.adj),ctot=P2.adj,m=m,...)
  y[1:R,1:R,]<-ipf$mu
  return(c(ipf,list(y=y)))
}
