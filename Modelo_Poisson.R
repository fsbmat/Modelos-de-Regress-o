setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

op <- par(mfrow = c(1, 1))
#f é uma função densidade de probabilidade
#Fixamos valores para os parâmetros, plotamos o gráfico da função e visualisamos que f>0 para 
#todo x>0 real e integramos f de 0 a +Infinito.

lambda=10
u=runif(10)

#Gráfico da densidade da distribuição Poisson 

plot(dpois(seq(1,10, by =1), lambda = lambda), type ="h",xlab =
       "Valores inteiros", ylab = "Probabilidade", main = "Função
     massa de probabilidade")

# #Gráfico da distribuição acumulada da Poisson 
# 
# plot(ppois(seq(1,10, by =1), lambda = lambda),type ="h", xlab =
#         "Valores Inteiros", ylab = "Probabilidade", main = "Função de
#       probabilidade acumulada")

#integral de f

integrand <- function(x) {dpois(x=seq(1,100000, by =1), lambda = lambda)}
sum(integrand(x))

#Note que a soma tende a 1

#Agora vamos simular valores de uma função Poisson

N=100
m=matrix(ncol=1,nrow=N)

for (i in 1:N){
  
lambda=2
n=100
x=rpois(n,lambda)

#hist(x,prob=T,ylim=c(0,1),xlim=c(0,20),breaks = 8, col = "grey", border = "black", main="Dados simulados com distribuição Poisson")
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

#Cálculo das médias de cada coluna da matriz de parâmetros m
mest=colMeans(m)

#Cálculo do desvio padrão de cada coluna da matriz de parâmetros m
dest=apply(m,2,sd)

#Cálculo do erro quadratíco médio de cada coluna da matriz de parâmetros m em relação ao verdadeiro valor do parâmetro
eqm=function(x,theta){ 
  N=length(x)
  sqrt(sum(((x-theta)^2))/N)}

#Erro quadrático médio estimado de cada um dos parâmetros 
eqmest=c(eqm(x=m[,1],theta=lambda))
eqmest

#Tabela com os verdadeiros valores dos parâmetros e com a média 
#desvio-padrão e erro quadrático médio dos parâmetros estimados
tab=data.frame(truevalue=c(lambda),mean=mest,sd=dest,eqm=eqmest)
tab

#Usando a Função Escore U(lambda)=-n+sum(x)/lambda. Como usei -log(verossimilhança), 
#multiplicamos U(lambda) por -1

UPois <- function(lambda, amostra){
  return(with(amostra, n - soma/lambda))
}

am <- list(n=length(x), soma=sum(x))
#Para obter a estimativa utilizamos inicialmente a função uniroot() que implementa
#um algoritmo para encontrar a raiz de uma equação.
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

#Função para simular variáveis aleatórias de uma modelo de regressão de Poisson.

simula.poisson <- function(formula, beta) {
  X <- model.matrix(formula)
  lambda <- exp(X %*% beta)
  y <- rpois(nrow(X), lambda = lambda)
#Tabela com os valores de y simulados e de x
    return(data.frame(y = y, X))
}



set.seed(123)
#Gerando valores para as covariáveis
cov <- seq(0, 5, length=10)
beta0=2
beta1=3
#Matriz dos dados simulados de y e com as covariáveis x, o vetor c() determina a quantidade de parâmetros a serem estimados
dados10 <- simula.poisson(~cov, c(beta0,beta1))

## Função escore
escore <- function(par, formula, dados){
  mf <- model.frame(formula, dados)
  X <- model.matrix(formula, data=mf)
  esco <- crossprod(model.response(mf) - exp(X %*% par), X)
  return(drop(esco))
}

## Hessiana ("naïve")
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

#Função de Newton Raphson

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

#Newton Raphson para estimação dos parâmetros Beta

#Chutes iniciais

beta0=beta0+0.5
beta1=beta1+0.5

(beta10 <- NewtonRaphson(initial=c(beta0,beta1), escore=escore,
                         hessiano=hessiano, max.iter=1000,
                         formula=y~cov, dados=dados10))


