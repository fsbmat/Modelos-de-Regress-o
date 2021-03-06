<<<<<<< HEAD
setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")
require(sos)


op <- par(mfrow = c(1, 1))

#X possui distribui��o exp(lambda)

#Gr�fico da densidade de X 

plot(function(x) (lambda*exp(-lambda*x)), 0,4)

#integral de f

integrand <- function(x) {lambda*exp(-lambda*x)}
integrate(integrand, lower = 0, upper = Inf)

#Note que a integral vale 1

#Agora vamos simular valores de uma vari�vel aleat�ria exponencial lambda

#M�todo de Monte Carlo

#N � o n�mero de amostras
N=1000
#Matriz onde vou armazenar todos os valores estimados ap�s cada la�o de Monte Carlo
m=matrix(ncol=1,nrow=N)

for (i in 1:N){
  
#Chute inicial do valor do par�metro 
  lambda=2
  #n � o tamanho de cada amostra
  n=100
  #Vamos gerar n valores da distribui��o uniforme(0,1)
  u=runif(n)
  #Fazemos X=F^(-1)(U) (F(x)=1-exp(-lambda*x)), assim X possui distribui��o exponencial lambda
  x=(-1/lambda)*log(1-u)
  #Histograma dos valores gerados
  #hist(x,prob=T,ylim=c(0,1),xlim=c(0,20),breaks = 8, col = "grey", border = "black", main="Dados simulados com distribui��o exponencial")
  #Curva te�rica da distribui��o exponencial(lambda)
  #curve(expr = (lambda*exp(-lambda*x)), from = 0.0001, to = 10,ylab = NULL,add=T)
  #Curva que acompanha os valores gerados
  #lines(density(x),col="blue")
  
#Defini��o da fun��o log de verossimilhan�a  
log.verossimilhanca <- function(par,x){ 
    lambda= par[1] 
    n=length(x)
    lv=(n*log(lambda)-lambda*sum(x)) 
    lv  
}
#Fun��o optim faz a estima��o do par�metro a partir da verossimilhan�a 
  theta=optim(par=2,x=x, fn=log.verossimilhanca, method="Nelder-Mead",hessian = F)
  #Fazendo a armazenagem de cada valor estimado em cada linha da matriz m
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

#Fun��o para simular vari�veis aleat�rias de uma modelo de regress�o exponencial.

simula.exponencial <- function(formula, beta) {
  X <- model.matrix(formula)
  lambda <- exp(-X %*% beta)
  y <- rexp(nrow(X), lambda)
  return(data.frame(y = y, X))
}

set.seed(123)
n=10
cov <- seq(0, 5, length=n)
dados1 <- simula.exponencial(~cov, c(2,0.5))
dados1

## Fun��o escore
escore <- function(par, formula, dados){
  mf <- model.frame(formula, dados)
  X <- model.matrix(formula, data=mf)
  esco <- crossprod((n*(t(X)*par)^(-1)-model.response(mf)), t(X))
  return(drop(esco))
  }

## Hessiana ("na�ve")
hessiano <- function(par, formula, dados){
  X <- model.matrix(formula, data=dados)
  mat <- matrix(0, nrow(X), nrow(X))
  diag(mat) <- -exp(X%*%par)
  H <- n/par^2
  return(H)
}


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


(beta10 <- NewtonRaphson(initial=c(0,0), escore=escore,
                         hessiano=hessiano, max.iter=1000,
                         formula=y~cov, dados=dados10))


=======
setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

op <- par(mfrow = c(1, 1))

#X possui distribui��o exp(lambda)

#Gr�fico da densidade de X 

plot(function(x) (lambda*exp(-lambda*x)), 0,4)

#integral de f

integrand <- function(x) {lambda*exp(-lambda*x)}
integrate(integrand, lower = 0, upper = Inf)

#Note que a integral vale 1

#Agora vamos simular valores de uma vari�vel aleat�ria exponencial lambda

#M�todo de Monte Carlo

#N � o n�mero de amostras
N=1000
#Matriz onde vou armazenar todos os valores estimados ap�s cada la�o de Monte Carlo
m=matrix(ncol=1,nrow=N)

for (i in 1:N){
  
#Chute inicial do valor do par�metro 
  lambda=2
  #n � o tamanho de cada amostra
  n=100
  #Vamos gerar n valores da distribui��o uniforme(0,1)
  u=runif(n)
  #Fazemos X=F^(-1)(U) (F(x)=1-exp(-lambda*x)), assim X possui distribui��o exponencial lambda
  x=(-1/lambda)*log(1-u)
  #Histograma dos valores gerados
  #hist(x,prob=T,ylim=c(0,1),xlim=c(0,20),breaks = 8, col = "grey", border = "black", main="Dados simulados com distribui��o exponencial")
  #Curva te�rica da distribui��o exponencial(lambda)
  #curve(expr = (lambda*exp(-lambda*x)), from = 0.0001, to = 10,ylab = NULL,add=T)
  #Curva que acompanha os valores gerados
  #lines(density(x),col="blue")
  
#Defini��o da fun��o log de verossimilhan�a  
log.verossimilhanca <- function(par,x){ 
    lambda= par[1] 
    n=length(x)
    lv=(n*log(lambda)-lambda*sum(x)) 
    lv  
}
#Fun��o optim faz a estima��o do par�metro a partir da verossimilhan�a 
  theta=optim(par=2,x=x, fn=log.verossimilhanca, method="Nelder-Mead",hessian = F)
  #Fazendo a armazenagem de cada valor estimado em cada linha da matriz m
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

#Fun��o para simular vari�veis aleat�rias de uma modelo de regress�o exponencial.

simula.exponencial <- function(formula, beta) {
  X <- model.matrix(formula)
  lambda <- exp(-X %*% beta)
  y <- rexp(nrow(X), lambda)
  return(data.frame(y = y, X))
}

set.seed(123)
n=10
cov <- seq(0, 5, length=n)
dados1 <- simula.exponencial(~cov, c(2,0.5))
dados1

## Fun��o escore
escore <- function(par, formula, dados){
  mf <- model.frame(formula, dados)
  X <- model.matrix(formula, data=mf)
  esco <- crossprod((n*(t(X)*par)^(-1)-model.response(mf)), t(X))
  return(drop(esco))
  }

## Hessiana ("na�ve")
hessiano <- function(par, formula, dados){
  X <- model.matrix(formula, data=dados)
  mat <- matrix(0, nrow(X), nrow(X))
  diag(mat) <- -exp(X%*%par)
  H <- n/par^2
  return(H)
}


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


(beta10 <- NewtonRaphson(initial=c(0,0), escore=escore,
                         hessiano=hessiano, max.iter=1000,
                         formula=y~cov, dados=dados10))


>>>>>>> origin/master
