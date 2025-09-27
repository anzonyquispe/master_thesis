# Dantzig for regression using quantreg

###################
# simulation design
###################

rm(list = ls())
setwd("~/Desktop/research/rrr_dantzig")
library("quantreg")

set.seed(1)
n <- 100
p <- 100
cf.true<- c(1, 1, 1, rep(0,p-2))

X <- matrix(rnorm(p*n), n, p)
X<- cbind(rep(1,p), X) # has a constant
p<- p+1


sigma = 1
y <- apply(X[, 1:3], 1, sum) + rnorm(n)*sigma 

G <- t(X) %*% X/n 
M <- c(as.matrix(t(X) %*% y))/n

get_mse<-function(y1,y2) {
  mse=mean((y1-y2)^2)
  return(mse)
}

two.norm <- function(x){
  return(sqrt(x %*% x))
} 


######################################################################################################################

# regression formulation

#######################################
# validating use of theoretical lambda
#######################################

DSelector <- function(X, y, sigma = 1, lambda = 3.5, sparse = TRUE)
{
  n <- nrow(X)
  p <- ncol(X)
  K <- lambda * sigma
  A <- t(X) %*% X
  R <- rbind(A, -A)
  a <- c(as.matrix(t(X) %*% y))
  r <- c(a - K, -a - K)
  zp <- rep(0, p)
  if(sparse){
    Ip <- as(p, "matrix.diag.csr")
    R <- as.matrix.csr(R)
    f <- rq.fit.sfnc(Ip, zp, R = R, r = r)
  }
  else{
    Ip <- diag(p)
    f <- rq.fit.fnc(Ip, zp, R = R, r = r)
  }
  return(f)
}


DSel <- function(X, y, D="Id", c=1, sparse = TRUE) { 
  n <- nrow(X) 
  p <- ncol(X) 
  if (is.matrix(D) ==0)  D= diag(rep(1,p))
  invD<- diag(1/diag(D))
  G <- invD%*%t(X) %*% X/n 
  R <- rbind(G, -G) 
  lambda<- (c)*qnorm(1-.1/(2*p))/sqrt(n)
  a <- invD%*%c(as.matrix(t(X) %*% y))/n
  r <- c(a - lambda, -a - lambda) 
  zp <- rep(0, p) 
  if(sparse){ Ip <- as(p, "matrix.diag.csr") 
  R <- as.matrix.csr(R) 
  f <- rq.fit.sfnc(Ip, zp, R = R, r = r) } 
  else{ Ip <- diag(p) 
  f <- rq.fit.fnc(Ip, zp, R = R, r = r) } 
  return(f) } 

cf.Koenker <- DSelector(X, y)$coef
cf.1 <- DSel(X, y,  D=diag(rep(sigma, p)), c= 1/2)$coef

MSE.Koenker = get_mse(cf.Koenker,cf.true) 
MSE.1 = get_mse(cf.1,cf.true) 

print(cbind(cf.Koenker, cf.1))            
print(cbind(MSE.Koenker, MSE.1))


######################################################################################################################

# RMD formulation, penalty on constant

###################################
# validating use of RMD formulation
###################################

D= diag(rep(1,p)) # can iterate to update it
c=0.5
alpha=0.1
lambda=c*qnorm(1-alpha/(2*p))/sqrt(n) 

RMD <- function(G, M, D, lambda, sparse = TRUE) { 
  p <- ncol(G) 
  invD<- diag(1/diag(D))
  Gt <- invD%*%G
  Mt <- invD%*%M
  R <- rbind(Gt, -Gt) 
  r <- c(Mt - lambda, -Mt - lambda) 
  zp <- rep(0, p) 
  if(sparse){ Ip <- as(p, "matrix.diag.csr") 
  R <- as.matrix.csr(R) 
  f <- rq.fit.sfnc(Ip, zp, R = R, r = r) } 
  else{ Ip <- diag(p) 
  f <- rq.fit.fnc(Ip, zp, R = R, r = r) } 
  return(f) } 

cf.2 <- RMD(G, M,  D, lambda)$coef
MSE.2 = get_mse(cf.2,cf.true) 

print(cbind(cf.Koenker, cf.2))            
print(cbind(MSE.Koenker, MSE.2))

######################################
# validating use of my RMD formulation
######################################

D= rep(1,p) # can iterate to update it

RMD_dantzig <- function(M, G, D, lambda=0, sparse = TRUE) {
  
  p <- ncol(G)
  zp <- rep(0, p)

  A <- solve(diag(D),G)
  R <- rbind(A, -A)
  
  a <- solve(diag(D),M)
  r <- c(a - lambda, -a - lambda)
  
  if(sparse) {
    Ip <- as(p, "matrix.diag.csr")
    R <- as.matrix.csr(R)
    f <- rq.fit.sfnc(Ip, zp, R = R, r = r)
  } else {
    Ip <- diag(p)
    f <- rq.fit.fnc(Ip, zp, R = R, r = r)
  }
  
  return(f)
}

cf.3 <- RMD_dantzig(M,G,D,lambda)$coef
MSE.3 = get_mse(cf.3,cf.true) 

print(cbind(cf.Koenker, cf.3))            
print(cbind(MSE.Koenker, MSE.3))

###########################
# validating use of scaling
###########################

get_D <- function(Y,X,rho_hat){
  n=nrow(X)
  p=ncol(X)
  
  df=matrix(0,p,n)
  for (i in 1:n){
    df[,i]=X[i,]*as.vector(rho_hat %*% X[i,])-Y[i]*X[i,]
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

# low-dimensional moments
p0=ceiling(p/2)
X0=X[,1:p0]
G0 <- t(X0) %*% X0/n 
M0 <- c(as.matrix(t(X0) %*% y))/n

# initial estimate
beta_hat0=solve(G0,M0)
beta_hat=c(beta_hat0,rep(0,p-p0))
beta_hat_old=beta_hat+0

# adaptive weight
D_hat_beta=get_D(y,X,beta_hat_old)
beta_hat=RMD_dantzig(M, G, D_hat_beta, lambda)$coefficients

cf.4a=beta_hat
MSE.4a = get_mse(cf.4a,cf.true) 

cf.4b=diag(D_hat_beta) %*% beta_hat # scaling by D as in pseudo-code
MSE.4b = get_mse(cf.4b,cf.true) 

print(cbind(cf.Koenker, cf.4a,cf.4b))            
print(cbind(MSE.Koenker, MSE.4a,MSE.4b)) #4a > Koenker > 4b

#############################
# validating use of iteration
#############################

# initial estimate
beta_hat0=solve(G0,M0)
beta_hat=c(beta_hat0,rep(0,p-p0))

tol=1e-6
diff=1

while(diff>tol){
  
  # previous values
  beta_hat_old=beta_hat+0
  
  # normalization
  D_hat_beta=get_D(y,X,beta_hat_old)
  
  # RMD estimate
  beta_hat=RMD_dantzig(M, G, D_hat_beta, lambda)$coefficients
  
  # difference
  diff=two.norm(beta_hat-beta_hat_old)
  
}

cf.5a=beta_hat
MSE.5a = get_mse(cf.5a,cf.true) 

cf.5b=diag(D_hat_beta) %*% beta_hat # scaling by D as in pseudo-code
MSE.5b = get_mse(cf.5b,cf.true) 

print(cbind(cf.Koenker, cf.5a,cf.5b))            
print(cbind(MSE.Koenker, MSE.5a,MSE.5b)) #4a > 5a > Koenker > 5b > 4b. 
# seems like we do not want to scale by D as final step. effect of iteration is ambiguous

################
# report results
################

print(round(cbind(MSE.Koenker,MSE.1,MSE.2,MSE.3,MSE.4b,MSE.4a,MSE.5b,MSE.5a),5))

######################################################################################################################

# RMD formulation, no penalty on constant

###################################
# validating use of RMD formulation
###################################

D= diag(rep(1,p)) # can iterate to update it
c=0.5
alpha=0.1
lambda=c*qnorm(1-alpha/(2*p))/sqrt(n) 
l=0.1

RMD2 <- function(G, M, D, lambda, sparse = TRUE) { 
  p <- ncol(G) 
  invD<- diag(1/diag(D))
  Gt <- invD%*%G
  Mt <- invD%*%M
  R <- rbind(Gt, -Gt) 
  L<- c(l, rep(1,p-1))  
  r <- c(Mt - lambda*L, -Mt - lambda*L) 
  zp <- rep(0, p) 
  if(sparse){ 
    Ip <- as(p, "matrix.diag.csr") 
    R <- as.matrix.csr(R) 
    f <- rq.fit.sfnc(Ip, zp, R = R, r = r) 
  }else{ 
    Ip <- diag(p) 
    f <- rq.fit.fnc(Ip, zp, R = R, r = r) 
  } 
  return(f) 
} 

cf.6 <- RMD2(G, M,  D, lambda)$coef
MSE.6 = get_mse(cf.6,cf.true) 

print(cbind(cf.Koenker, cf.6))            
print(cbind(MSE.Koenker, MSE.6))


######################################
# validating use of my RMD formulation
######################################

D= rep(1,p) # can iterate to update it

RMD_dantzig2 <- function(M, G, D, lambda=0, sparse = TRUE) {
  
  p <- ncol(G)
  zp <- rep(0, p)
  L<- c(l, rep(1,p-1))  
  
  A <- solve(diag(D),G)
  R <- rbind(A, -A)
  
  a <- solve(diag(D),M)
  r <- c(a - lambda*L, -a - lambda*L)
  
  if(sparse) {
    Ip <- as(p, "matrix.diag.csr")
    R <- as.matrix.csr(R)
    f <- rq.fit.sfnc(Ip, zp, R = R, r = r)
  } else {
    Ip <- diag(p)
    f <- rq.fit.fnc(Ip, zp, R = R, r = r)
  }
  
  return(f)
}

cf.7 <- RMD_dantzig2(M,G,D,lambda)$coef
MSE.7 = get_mse(cf.7,cf.true) 

print(cbind(cf.Koenker, cf.7))            
print(cbind(MSE.Koenker, MSE.7))

###########################
# validating use of scaling
###########################

get_D <- function(Y,X,rho_hat){
  n=nrow(X)
  p=ncol(X)
  
  df=matrix(0,p,n)
  for (i in 1:n){
    df[,i]=X[i,]*as.vector(rho_hat %*% X[i,])-Y[i]*X[i,]
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

# low-dimensional moments
p0=ceiling(p/2)
X0=X[,1:p0]
G0 <- t(X0) %*% X0/n 
M0 <- c(as.matrix(t(X0) %*% y))/n

# initial estimate
beta_hat0=solve(G0,M0)
beta_hat=c(beta_hat0,rep(0,p-p0))
beta_hat_old=beta_hat+0

# adaptive weight
D_hat_beta=get_D(y,X,beta_hat_old)
beta_hat=RMD_dantzig2(M, G, D_hat_beta, lambda)$coefficients

cf.8a=beta_hat
MSE.8a = get_mse(cf.8a,cf.true) 

cf.8b=diag(D_hat_beta) %*% beta_hat # scaling by D as in pseudo-code
MSE.8b = get_mse(cf.8b,cf.true) 

print(cbind(cf.Koenker, cf.8a,cf.8b))            
print(cbind(MSE.Koenker, MSE.8a,MSE.8b))

#############################
# validating use of iteration
#############################

# initial estimate
beta_hat0=solve(G0,M0)
beta_hat=c(beta_hat0,rep(0,p-p0))

tol=1e-6
diff=1

while(diff>tol){
  
  # previous values
  beta_hat_old=beta_hat+0
  
  # normalization
  D_hat_beta=get_D(y,X,beta_hat_old)
  
  # RMD estimate
  beta_hat=RMD_dantzig2(M, G, D_hat_beta, lambda)$coefficients
  
  # difference
  diff=two.norm(beta_hat-beta_hat_old)
  
}

cf.9a=beta_hat
MSE.9a = get_mse(cf.9a,cf.true) 

cf.9b=diag(D_hat_beta) %*% beta_hat # scaling by D as in pseudo-code
MSE.9b = get_mse(cf.9b,cf.true) 

print(cbind(cf.Koenker, cf.9a,cf.9b))            
print(cbind(MSE.Koenker, MSE.9a,MSE.9b))

################
# report results
################

print(round(cbind(MSE.6,MSE.7,MSE.8b,MSE.8a,MSE.9b,MSE.9a),5))



