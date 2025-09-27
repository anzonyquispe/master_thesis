two.norm <- function(x){
  return(sqrt(x %*% x))
} 

one.norm<-function(x){
  return(sum(x%*%sign(x)))
}

one.norm.grad<-function(x){
  return(sign(x))
}

m<-function(y,x,dx,beta){ #all data arguments to make interchangeable with m2
  return(dx %*% beta)
}

m2<-function(y,x,dx,beta){
  return(y*as.vector(x %*% beta))
}

psi_tilde<-function(y,x,dx,m,rho,beta){
  return(m(y,x,dx,beta)+(x %*% rho)*(y-(x %*% beta)))
}

psi_tilde_bias<-function(y,x,dx,m,rho,beta){
  return(m(y,x,dx,beta))
}

get_MNG<-function(Y,X,dX){
  
  p=ncol(X)
  n=nrow(X)

  N=Y*X
  
  M_hat=colMeans(dX) #since m(w,b)=dx. do we want equal weighting of individuals? yes
  N_hat=colMeans(N)
  G_hat=t(X)%*%X/n
  
  return(list(M_hat,N_hat,G_hat))
}