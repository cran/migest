rc9 <-
function(x, param=rc9.fund, scaled=TRUE){
  if(!is.list(param))
    stop("param must be a list")
  if(sum(is.na(match(names(param),names(rc9.fund))))!=0)
    stop("param must be a list with correct names, see for example rc9.fund")
  m<-param$a1*exp(-param$alpha1*x)+param$a2*exp(-param$alpha2*(x-param$mu2)-exp(-param$lambda*(x-param$mu2)))+param$c
  if(scaled==TRUE) m<-m/sum(m)
  m
}
