setwd("~/GitHub/Modelos-de-Regress-o")
rm(list=ls())
cat("\014")

op <- par(mfrow = c(1, 1))

#X possui distribuição exp(lambda)

#Gráfico da densidade de X 

plot(function(x) (lambda*exp(-lambda*x)), 0,4)

#integral de f

integrand <- function(x) {lambda*exp(-lambda*x)}
integrate(integrand, lower = 0, upper = Inf)

#Note que a integral vale 1

#Agora vamos simular valores de uma variável aleatória exponencial lambda

#Método de Monte Carlo

#N é o número de amostras
N=1000
#Matriz onde vou armazenar todos os valores estimados após cada laço de Monte Carlo
m=matrix(ncol=1,nrow=N)

for (i in 1:N){
  
#Chute inicial do valor do parâmetro 
  lambda=2
  #n é o tamanho de cada amostra
  n=100
  #Vamos gerar n valores da distribuição uniforme(0,1)
  u=runif(n)
  #Fazemos X=F^(-1)(U) (F(x)=1-exp(-lambda*x)), assim X possui distribuição exponencial lambda
  x=(-1/lambda)*log(1-u)
  #Histograma dos valores gerados
  #hist(x,prob=T,ylim=c(0,1),xlim=c(0,20),breaks = 8, col = "grey", border = "black", main="Dados simulados com distribuição exponencial")
  #Curva teórica da distribuição exponencial(lambda)
  #curve(expr = (lambda*exp(-lambda*x)), from = 0.0001, to = 10,ylab = NULL,add=T)
  #Curva que acompanha os valores gerados
  #lines(density(x),col="blue")
  
#Definição da função log de verossimilhança  
log.verossimilhanca <- function(par,x){ 
    lambda= par[1] 
    n=length(x)
    lv=(n*log(lambda)-lambda*sum(x)) 
    lv  
}
#Função optim faz a estimação do parâmetro a partir da verossimilhança 
  theta=optim(par=2,x=x, fn=log.verossimilhanca, method="Nelder-Mead",hessian = F)
  #Fazendo a armazenagem de cada valor estimado em cada linha da matriz m
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

#Função para simular variáveis aleatórias de uma modelo de regressão exponencial.

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

## Função escore
escore <- function(par, formula, dados){
  mf <- model.frame(formula, dados)
  X <- model.matrix(formula, data=mf)
  esco <- crossprod((n*(t(X)*par)^(-1)-model.response(mf)), t(X))
  return(drop(esco))
  }

## Hessiana ("naïve")
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


