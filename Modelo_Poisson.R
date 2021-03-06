<<<<<<< HEAD
setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

op <- par(mfrow = c(1, 1))
#f � uma fun��o densidade de probabilidade
#Fixamos valores para os par�metros, plotamos o gr�fico da fun��o e visualisamos que f>0 para 
#todo x>0 real e integramos f de 0 a +Infinito.

lambda=10
u=runif(10)

#Gr�fico da densidade da distribui��o Poisson 

plot(dpois(seq(1,10, by =1), lambda = lambda), type ="h",xlab =
       "Valores inteiros", ylab = "Probabilidade", main = "Fun��o
     massa de probabilidade")

# #Gr�fico da distribui��o acumulada da Poisson 
# 
# plot(ppois(seq(1,10, by =1), lambda = lambda),type ="h", xlab =
#         "Valores Inteiros", ylab = "Probabilidade", main = "Fun��o de
#       probabilidade acumulada")

#integral de f

integrand <- function(x) {dpois(x=seq(1,100000, by =1), lambda = lambda)}
sum(integrand(x))

#Note que a soma tende a 1

#Agora vamos simular valores de uma fun��o Poisson

N=100
m=matrix(ncol=1,nrow=N)

for (i in 1:N){
  
lambda=2
n=100
x=rpois(n,lambda)

#hist(x,prob=T,ylim=c(0,1),xlim=c(0,20),breaks = 8, col = "grey", border = "black", main="Dados simulados com distribui��o Poisson")
#curve(expr=(ppois(x=seq(1,10, by =1), lambda = lambda), from = 0, to = 10,ylab = NULL,add=T)
#lines(density(x),col="blue")
  
  
log.verossimilhanca <- function(par, x){
   lambda= par[1]
   n=length(x)
   lv <- -sum(dpois(x, lambda = lambda, log = TRUE))
   return(lv)
}
theta=optim(par=1,x=x, fn=log.verossimilhanca, method="Nelder-Mead",hessian = F) 
m[i,]=theta$par
}
#m

#C�lculo das m�dias de cada coluna da matriz de par�metros m
mest=colMeans(m)

#C�lculo do desvio padr�o de cada coluna da matriz de par�metros m
dest=apply(m,2,sd)

#C�lculo do erro quadrat�co m�dio de cada coluna da matriz de par�metros m em rela��o ao verdadeiro valor do par�metro
eqm=function(x,theta){ 
  N=length(x)
  sqrt(sum(((x-theta)^2))/N)}

#Erro quadr�tico m�dio estimado de cada um dos par�metros 
eqmest=c(eqm(x=m[,1],theta=lambda))
eqmest

#Tabela com os verdadeiros valores dos par�metros e com a m�dia 
#desvio-padr�o e erro quadr�tico m�dio dos par�metros estimados
tab=data.frame(truevalue=c(lambda),mean=mest,sd=dest,eqm=eqmest)
tab

#Usando a Fun��o Escore U(lambda)=-n+sum(x)/lambda. Como usei -log(verossimilhan�a), 
#multiplicamos U(lambda) por -1

UPois <- function(lambda, amostra){
  return(with(amostra, n - soma/lambda))
}

am <- list(n=length(x), soma=sum(x))
#Para obter a estimativa utilizamos inicialmente a fun��o uniroot() que implementa
#um algoritmo para encontrar a raiz de uma equa��o.
uniroot(UPois, interval=range(x), amostra=am)$root

#Usando Newton Raphson, precisamos encontrar a matriz Hessiana que denominei HPois

HPois <- function(lambda, amostra){
  return(amostra$soma/lambda^2)
}

maxit <- 100; lambdaNR <- lambda; iter <- 0; d <- 1
while(d > 1e-12 & iter <= maxit){
  lambdaNR.new <-
  lambdaNR - UPois(lambdaNR, am)/HPois(lambdaNR, am)
  d <- abs(lambdaNR - lambdaNR.new)
  lambdaNR <- lambdaNR.new ; iter <- iter + 1
}
c(lambdaNR, iter)

########################################################################

#Fun��o para simular vari�veis aleat�rias de uma modelo de regress�o de Poisson.

simula.poisson <- function(formula, beta) {
  X <- model.matrix(formula)
  lambda <- exp(X %*% beta)
  y <- rpois(nrow(X), lambda = lambda)
#Tabela com os valores de y simulados e de x
    return(data.frame(y = y, X))
}



set.seed(123)
#Gerando valores para as covari�veis
cov <- seq(0, 5, length=10)
beta0=2
beta1=3
#Matriz dos dados simulados de y e com as covari�veis x, o vetor c() determina a quantidade de par�metros a serem estimados
dados10 <- simula.poisson(~cov, c(beta0,beta1))

## Fun��o escore
escore <- function(par, formula, dados){
  mf <- model.frame(formula, dados)
  X <- model.matrix(formula, data=mf)
  esco <- crossprod(model.response(mf) - exp(X %*% par), X)
  return(drop(esco))
}

## Hessiana ("na�ve")
hessiano <- function(par, formula, dados){
  X <- model.matrix(formula, data=dados)
  mat <- matrix(0, nrow(X), nrow(X))
  diag(mat) <- -exp(X%*%par)
  H <- t(X) %*% mat %*% X
  return(H)
}

## Hessiana (equivalente a anterior)
hessiano <- function(par, formula, dados){
  X <- model.matrix(formula, data=dados)
  H <- crossprod(X * -exp(drop(X%*%par)), X)
  return(H)
}

#Fun��o de Newton Raphson

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

#Newton Raphson para estima��o dos par�metros Beta

#Chutes iniciais

beta0=beta0+0.5
beta1=beta1+0.5

(beta10 <- NewtonRaphson(initial=c(beta0,beta1), escore=escore,
                         hessiano=hessiano, max.iter=1000,
                         formula=y~cov, dados=dados10))


=======
setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

op <- par(mfrow = c(1, 1))
#f � uma fun��o densidade de probabilidade
#Fixamos valores para os par�metros, plotamos o gr�fico da fun��o e visualisamos que f>0 para 
#todo x>0 real e integramos f de 0 a +Infinito.

lambda=10
u=runif(10)

#Gr�fico da densidade da distribui��o Poisson 

plot(dpois(seq(1,10, by =1), lambda = lambda), type ="h",xlab =
       "Valores inteiros", ylab = "Probabilidade", main = "Fun��o
     massa de probabilidade")

#Gr�fico da distribui��o acumulada da Poisson Poisson 

plot(ppois(seq(1,10, by =1), lambda = lambda),type ="h", xlab =
        "Valores Inteiros", ylab = "Probabilidade", main = "Fun��o de
      probabilidade acumulada")

#integral de f

integrand <- function(x) {dpois(x=seq(1,100000, by =1), lambda = lambda)}
sum(integrand(x))

#Note que a soma tende a 1

#Agora vamos simular valores de uma fun��o Poisson

N=100
m=matrix(ncol=1,nrow=N)

for (i in 1:N){
  
  lambda=2
  n=100
  x=rpois(n,lambda)
  #hist(x,prob=T,ylim=c(0,1),xlim=c(0,20),breaks = 8, col = "grey", border = "black", main="Dados simulados com distribui��o Poisson")
  #curve(expr=(ppois(x=seq(1,10, by =1), lambda = lambda), from = 0, to = 10,ylab = NULL,add=T)
  #lines(density(x),col="blue")
  
  
 log.verossimilhanca <- function(par, x){
   lambda= par[1]
   n=length(x)
   lv <- -sum(dpois(x, lambda = lambda, log = TRUE))
   return(lv)
 }
  theta=optim(par=1,x=x, fn=log.verossimilhanca, method="Nelder-Mead",hessian = F) 
  m[i,]=theta$par
}
#m

#C�lculo das m�dias de cada coluna da matriz de par�metros m
mest=colMeans(m)

#C�lculo do desvio padr�o de cada coluna da matriz de par�metros m
dest=apply(m,2,sd)

#C�lculo do erro quadrat�co m�dio de cada coluna da matriz de par�metros m em rela��o ao verdadeiro valor do par�metro
eqm=function(x,theta){ 
  N=length(x)
  sqrt(sum(((x-theta)^2))/N)}

#Erro quadr�tico m�dio estimado de cada um dos par�metros 
eqmest=c(eqm(x=m[,1],theta=lambda))
eqmest

#Tabela com os verdadeiros valores dos par�metros e com a m�dia 
#desvio-padr�o e erro quadr�tico m�dio dos par�metros estimados
tab=data.frame(truevalue=c(lambda),mean=mest,sd=dest,eqm=eqmest)
tab

#Usando a Fun��o Escore U(lambda)=-n+sum(x)/lambda. Como usei -log(verossimilhan�a), 
#multiplicamos U(lambda) por -1

UPois <- function(lambda, amostra){
  return(with(amostra, n - soma/lambda))
}

am <- list(n=length(x), soma=sum(x))
#Para obter a estimativa utilizamos inicialmente a fun��o uniroot() que implementa
#um algoritmo para encontrar a raiz de uma equa��o.
uniroot(UPois, interval=range(x), amostra=am)$root

#Usando Newton Raphson, precisamos encontrar a matriz Hessiana que denominei HPois

HPois <- function(lambda, amostra){
  return(amostra$soma/lambda^2)
}

maxit <- 100; lambdaNR <- lambda; iter <- 0; d <- 1
while(d > 1e-12 & iter <= maxit){
  lambdaNR.new <-
  lambdaNR - UPois(lambdaNR, am)/HPois(lambdaNR, am)
  d <- abs(lambdaNR - lambdaNR.new)
  lambdaNR <- lambdaNR.new ; iter <- iter + 1
}
c(lambdaNR, iter)

########################################################################

#Fun��o para simular vari�veis aleat�rias de uma modelo de regress�o de Poisson.

simula.poisson <- function(formula, beta) {
  X <- model.matrix(formula)
  lambda <- exp(X %*% beta)
  y <- rpois(nrow(X), lambda = lambda)
#Tabela com os valores de y simulados e de x
    return(data.frame(y = y, X))
}



set.seed(123)
#Gerando valores para as covari�veis
cov <- seq(0, 5, length=10)
beta0=2
beta1=3
#Matriz dos dados simulados de y e com as covari�veis x, o vetor c() determina a quantidade de par�metros a serem estimados
dados10 <- simula.poisson(~cov, c(beta0,beta1))

## Fun��o escore
escore <- function(par, formula, dados){
  mf <- model.frame(formula, dados)
  X <- model.matrix(formula, data=mf)
  esco <- crossprod(model.response(mf) - exp(X %*% par), X)
  return(drop(esco))
}

## Hessiana ("na�ve")
hessiano <- function(par, formula, dados){
  X <- model.matrix(formula, data=dados)
  mat <- matrix(0, nrow(X), nrow(X))
  diag(mat) <- -exp(X%*%par)
  H <- t(X) %*% mat %*% X
  return(H)
}

## Hessiana (equivalente a anterior)
hessiano <- function(par, formula, dados){
  X <- model.matrix(formula, data=dados)
  H <- crossprod(X * -exp(drop(X%*%par)), X)
  return(H)
}

#Fun��o de Newton Raphson

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

#Newton Raphson para estima��o dos par�metros Beta

#Chutes iniciais

beta0=beta0+0.5
beta1=beta1+0.5

(beta10 <- NewtonRaphson(initial=c(beta0,beta1), escore=escore,
                         hessiano=hessiano, max.iter=1000,
                         formula=y~cov, dados=dados10))


>>>>>>> origin/master
