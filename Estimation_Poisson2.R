N=1000
m=matrix(ncol=2,nrow=N)
m1=matrix(ncol=2,nrow=N)
m2=matrix(ncol=2,nrow=N)
for (i in 1:N){
#set.seed(123)  
beta0=2
beta1=0.5
n=5
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
  
hess<-function(beta){
  lambda=exp(X%*%beta[1:p])
  hess11=-sum(lambda)
  hess12=-sum(t(X1)%*%lambda)
  hess21=-sum(t(X1)%*%lambda)
  hess22=-sum(t(X1)%*%lambda%*%t(X1))
  matrix(c(-hess11,-hess21, -hess12,-hess22),nrow=2,ncol=2)
}
  
  
const1 <- rep(1,n);
const <- cbind(const1);
X1=matrix(runif(n, 0, 1),nrow=n,ncol=1)
X <-matrix(c(const,X1),nrow=n,ncol=ncol(X1)+1)
p=ncol(X)
start=c(2,2)

#Estimation with Newton Raphson
NewtonRaphson <- function(initial, escore, hessiano, tol=0.0001,
                            max.iter, n.dim, print=FALSE, ...){
solucao <- initial
  for(i in 2:max.iter){
    solucao <- initial - solve(hessiano(initial, ...),
                             escore(initial, ...))
    tolera <- abs(solucao - initial)
    if(all(tolera < tol) == TRUE)break
    initial <- solucao
    if(print) print(initial)
 }
 return(initial)
}

beta10 <- NewtonRaphson(initial=c(2,0.5), escore=grad,hessiano=hess, max.iter=1000)
m2[i,]=beta10  

#Estimation with function optim
poisson_opt=optim(start,like_poisson,grad,method="BFGS",hessian = T) 
m[i,]=poisson_opt$par
#Estimation with function nlminb
poisson_nlm=nlminb(start, like_poisson,grad)
m1[i,]=poisson_nlm$par
}

#
#Calculating the average of each column of the array of parameters m
mest=colMeans(m)
mest1=colMeans(m1)
mest2=colMeans(m2)
#calculating the standard deviation of each column of the array of parameters m
dest=apply(m,2,sd)
dest1=apply(m1,2,sd)
dest2=apply(m2,2,sd)
#root mean square error in the calculation of each column of the array of parameters m in relation to the true value of the parameter
eqm=function(x,poisson_opt){ 
k=length(x)
sqrt(sum(((x-poisson_opt)^2))/k)}

eqm1=function(x,poisson_nlm){ 
k=length(x)
sqrt(sum(((x-poisson_nlm)^2))/k)}

eqm2=function(x,beta10){ 
  k=length(x)
  sqrt(sum(((x-beta10)^2))/k)}

#Estimated mean squared error of each parameter 
beta0=2
beta1=0.5
lambda=c(beta0,beta1)
eqmest=c(eqm(x=m[,1],poisson_opt=lambda[1]),
         eqm(x=m[,2],poisson_opt=lambda[2]))

#Estimated mean squared error of each parameter
eqmest1=c(eqm1(x=m1[,1],poisson_nlm=lambda[1]),
          eqm1(x=m1[,2],poisson_nlm=lambda[2]))

#Estimated mean squared error of each parameter
eqmest2=c(eqm2(x=m2[,1],beta10=lambda[1]),
          eqm2(x=m2[,2],beta10=lambda[2]))


# Table with the true values of the parameters and the average
# Standard deviation and mean square error of the estimated parameters
tab=data.frame(truevalue=lambda,mean=mest,sd=dest,eqm=eqmest)
tab1=data.frame(truevalue=lambda,mean=mest1,sd=dest1,eqm=eqmest1)
tab2=data.frame(truevalue=lambda,mean=mest2,sd=dest2,eqm=eqmest2)
tab
tab1
tab2
