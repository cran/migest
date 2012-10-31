ffs <-
function(P1,P2,d,b,m=NULL,...){
  R<-nrow(P1)
  #set up offset
  if(is.null(m)){
    m<-array(1,c(dim(P1),dim(P1)[1]))
  }
  if(is.null(dimnames(m))){
    dimnames(m)<-list(orig=dimnames(P1)[[1]],dest=dimnames(P2)[[1]],pob=dimnames(P1)[[2]])
  }
  #set up migration matrix with births/deaths and others rows and columns
  y<-m
  y<-array(0,dim(m)+c(2,2,0))
  dimnames(y)<-list(o=c(dimnames(m)[[1]],"b","O"),d=c(dimnames(m)[[1]],"d","O"),b=dimnames(m)[[3]])
  y[,,]<-0
  
  #step 1-2a take off deaths
  d.mat<-t(ipf2(rt=d,m=t(P1))$mu)
  #d.mat<-round(d.mat)
  y[1:R,R+1,]<-t(d.mat)
  P1.adj<-P1-d.mat
  
  #step 1-2b take off births
  b.mat<-diag(b)
  y[R+1,1:R,]<-b.mat
  P2.adj<-P2-b.mat
  
  #step 3-4a take off moves in from external 
  dif<-rowSums(P1.adj) - rowSums(P2.adj)
  in.mat<-t(ipf2(ct=pmax(dif,0),m=t(P2.adj))$mu)
  P1.adj<-P1.adj-in.mat
  y[R+2,1:R,]<-t(in.mat)
  
  #step 3-4b take off moves out from external
  out.mat<-t(ipf2(ct=pmax(-dif,0),m=t(P1.adj))$mu)
  #out.mat<-round(out.mat)
  P2.adj<-P2.adj-out.mat
  y[1:R,R+2,]<-t(out.mat)
  
  #step 5 calculate
  ipf<-ipf3.qi(rt=t(P1.adj),ct=P2.adj,m=m,...)
  y[1:R,1:R,]<-ipf$mu
  return(c(ipf,list(y=y)))
}
