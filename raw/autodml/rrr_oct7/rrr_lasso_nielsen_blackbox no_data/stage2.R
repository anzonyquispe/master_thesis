L=5

rrr<-function(hh,Y,X,dX,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator,bias){
  
  N=length(unique(hh)) #how many hh
  sum_Ti=length(hh) # how many obs
  
  folds <- split(sample(N, N,replace=FALSE), as.factor(1:L)) # preserve hh in fold for clustering
  
  Psi_tilde=numeric(0)
  hh_new=numeric(0) #since sample splitting will re-order influences, with consequences for clustering of unbalanced panel
  
  for (l in 1:L){
    
    hh.l=folds[[l]] # which hh in this fold
    hh_idx.l=hh[hh %in% hh.l] # which hh in this fold, with repeats
    N.l=length(hh.l) # how many hh in this fold
    sum_Ti.l=length(hh_idx.l) # how many obs in this fold

    Y.l=Y[hh %in% hh.l]
    Y.nl=Y[!(hh %in% hh.l)]
    
    X.l=X[hh %in% hh.l,]
    X.nl=X[!(hh %in% hh.l),]
    
    dX.l=dX[hh %in% hh.l,]
    dX.nl=dX[!(hh %in% hh.l),]
    
    n.l=length(Y.l)
    n.nl=length(Y.nl)
    
    # get stage 1 (on nl)
    stage1_estimators<-get_stage1(Y.nl,X.nl,dX.nl,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator)
    rho_hat=stage1_estimators[[1]]
    beta_hat=stage1_estimators[[2]]
    
    # debugging
    print(paste0('fold: ',l))
    #print(paste0('beta_hat: '))
    #print(paste0(round(beta_hat,2)))
    #print(paste0('rho_hat: '))
    #print(paste0(round(rho_hat,2)))
    
    #get stage 2 (on l)
    Psi_tilde.l=rep(0,n.l)
    for (i in 1:n.l){
      if(bias){ #plug-in
        Psi_tilde.l[i]=psi_tilde_bias(Y.l[i],X.l[i,],dX.l[i,],m,rho_hat,beta_hat) # without subtracting theta_hat
      }else{ #DML
        Psi_tilde.l[i]=psi_tilde(Y.l[i],X.l[i,],dX.l[i,],m,rho_hat,beta_hat) # without subtracting theta_hat
      }
    }
    
    Psi_tilde=c(Psi_tilde,Psi_tilde.l)
    hh_new=c(hh_new,hh_idx.l) #hh indexing for Psi_tilde and hence Psi
    
  }
  
  # point estimate
  bar_Y=mean(Y)
  
  ape=mean(Psi_tilde)
  elasticity=ape/bar_Y-1
  
  # influences
  Psi=Psi_tilde-ape
  
  #####################
  # clustered variances
  #####################
  
  # vector with a component for each hh
  to_sum11<-rep(0,N)
  to_sum12<-rep(0,N)
  to_sum22<-rep(0,N)
  
  #for each hh j obtain sum of outer product
  for (j in 1:N){ 
    Psi.j=Psi[hh_new==j] #Psi for hh j
    Y.j=Y[hh==j] #Y for hh j
    Y.j_centered=Y.j-mean(Y) #centered Y for hh j
    
    to_sum11[j]=sum(Psi.j %o% Psi.j)
    to_sum12[j]=sum(Psi.j %o% Y.j)
    to_sum22[j]=sum(Y.j_centered %o% Y.j_centered)
  }
  
  # V
  v11=sum(to_sum11)/sum_Ti
  v12=sum(to_sum12)/sum_Ti  
  v21=v12
  v22=sum(to_sum22)/sum_Ti
  V=matrix(c(v11,v21,v12,v22),nrow=2)
  
  # H
  h1=bar_Y
  h2=-ape/(bar_Y^2)
  H=matrix(c(h1,h2),nrow=1)
  
  # delta method
  var=H %*% V %*% t(H)
  se=sqrt(var/sum_Ti)
  out<-c(sum_Ti,elasticity,se)
  
  return(out)
}
printer<-function(spec1){
  print(paste("   n:    ",spec1[1],"   elasticity:    ",round(spec1[2],2), "   SE:   ", round(spec1[3],2), sep=""))
}

for_tex<-function(spec1){
  print(paste(" & ",spec1[1]," & ",round(spec1[2],2), "   &   ", round(spec1[3],2), sep=""))
}