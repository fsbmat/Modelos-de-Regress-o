setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

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

#####################################################################################
####                             Modelo de Regressão                            #####
#####################################################################################

setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

N=100
m=matrix(ncol=1,nrow=N)
m1=matrix(ncol=1,nrow=N)
for (i in 1:N){

lambda=2
n=100
y=rpois(n,lambda)
remove(lambda)
#P(Yi=y)=(exp(-lambda)*lambda^y)/y!

like_poisson<- function(beta){
  lambda=exp(X%*%beta)
  loglik <- -sum(-lambda+y*log(lambda)-log(factorial(y)))
  return(loglik)
}

grad<-function(beta){
  lambda=exp(X%*%beta)
  gd=sum(-1+y/lambda)
  return(eval(-gd))
}

#####################################
#######Modelo Nulo beta=(beta0)######
#####################################

const <- rep(1,length(y));
X <- cbind(const);
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
tab=data.frame(truevalue=log(lambda),mean=mest,sd=dest,eqm=eqmest)
tab1=data.frame(truevalue=log(lambda),mean=mest1,sd=dest1,eqm=eqmest1)
tab
tab1

##########################################
#######Modelo com beta=(beta0,beta1)######
##########################################

setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

N=1000
m=matrix(ncol=2,nrow=N)
m1=matrix(ncol=2,nrow=N)
for (i in 1:N){
#set.seed(123)  
beta0=2
beta1=0.5
n=100
beta=matrix(c(beta0,beta1),nrow=2,ncol=1)
const1 <- rep(1,n);
const <- cbind(const1);
X1=matrix(runif(n, 0, 1),nrow=n,ncol=1)
X <-matrix(c(const,X1),nrow=n,ncol=ncol(X1)+1)
lambda=exp(X%*%beta)
y=rpois(n,lambda)
remove(beta0,beta1,X1,X)
#P(Yi=y)=(exp(-lambda)*lambda^y)/y!

like_poisson<- function(beta){
    lambda=exp(X%*%beta[1:p])
    loglik <- -sum(-lambda+y*log(lambda)-log(factorial(y)))
    return(loglik)
}

grad<-function(beta){
    lambda=exp(X%*%beta[1:p])
    grad1=sum(-lambda+y)
    grad2=sum(-t(X1)%*%lambda+y%*%X1)
    c(-grad1,-grad2)
}

const1 <- rep(1,n);
const <- cbind(const1);
X1=matrix(runif(n, 0, 1),nrow=n,ncol=1)
X <-matrix(c(const,X1),nrow=n,ncol=ncol(X1)+1)
p=ncol(X)
start=c(2,2)
  
#function optim
poisson_opt=optim(start,like_poisson,grad,method="BFGS",hessian = T) 
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
  k=length(x)
  sqrt(sum(((x-poisson_opt)^2))/k)}

eqm1=function(x,poisson_nlm){ 
  k=length(x)
  sqrt(sum(((x-poisson_nlm)^2))/k)}

#Estimated mean squared error of each parameter 
beta0=2
beta1=0.5
lambda=c(beta0,beta1)
eqmest=c(eqm(x=m[,1],poisson_opt=lambda[1]),
         eqm(x=m[,2],poisson_opt=lambda[2]))

#Estimated mean squared error of each parameter
eqmest1=c(eqm1(x=m1[,1],poisson_nlm=lambda[1]),
          eqm1(x=m1[,2],poisson_nlm=lambda[2]))


# Table with the true values of the parameters and the average
# Standard deviation and mean square error of the estimated parameters
tab=data.frame(truevalue=lambda,mean=mest,sd=dest,eqm=eqmest)
tab1=data.frame(truevalue=lambda,mean=mest1,sd=dest1,eqm=eqmest1)
tab
tab1

