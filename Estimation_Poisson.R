rm(list=ls())

N=100
m=matrix(ncol=1,nrow=N)
m1=matrix(ncol=1,nrow=N)
for (i in 1:N){
  
  lambda=2
  n=100
  y=rpois(n,lambda)
  remove(lambda)
  #P(Yi=y)=(exp(-lambda)*lambda^y)/y!
  
  like_poisson<- function(lambda){
    loglik <- -sum(-lambda+y*log(lambda)-log(factorial(y)))
    return(loglik)
  }
  
  grad<-function(lambda){
    gd=sum(-1+y/lambda)
    return(eval(-gd))
  }
  
  start=0.1
  
  #function optim
  poisson_opt=optim(start,like_poisson, method="BFGS",grad,hessian = F) 
  m[i,]=poisson_opt$par
  #function nlminb
  poisson_nlm=nlminb(start, like_poisson,grad)
  m1[i,]=poisson_nlm$par
}

#
#Calculating the average of each column of the array of parameters m
mest=colMeans(m)
mest1=colMeans(m1)

#calculating the standard deviation of each column of the array of parameters m
dest=apply(m,2,sd)
dest1=apply(m1,2,sd)
#root mean square error in the calculation of each column of the array of parameters m in relation to the true value of the parameter
eqm=function(x,poisson_opt){ 
  N=length(x)
  sqrt(sum(((x-poisson_opt)^2))/N)}

eqm1=function(x,poisson_nlm){ 
  N=length(x)
  sqrt(sum(((x-poisson_nlm)^2))/N)}

#Estimated mean squared error of each parameter 
lambda=2
eqmest=c(eqm(x=m[,1],poisson_opt=lambda))

#Estimated mean squared error of each parameter
eqmest1=c(eqm1(x=m1[,1],poisson_nlm=lambda))


# Table with the true values of the parameters and the average
# Standard deviation and mean square error of the estimated parameters
tab=data.frame(truevalue=lambda,mean=mest,sd=dest,eqm=eqmest)
tab1=data.frame(truevalue=lambda,mean=mest1,sd=dest1,eqm=eqmest1)
tab
tab1