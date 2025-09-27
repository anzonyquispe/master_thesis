L=5

rrr<-function(Y,T,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias,c_RR,c_CEF){
  
  n=nrow(X)
  folds <- split(sample(n, n,replace=FALSE), as.factor(1:L))
  
  Psi_tilde=numeric(0)
  Y_new=numeric(0)
  T_new=numeric(0)
  
  for (l in 1:L){
    
    Y.l=Y[folds[[l]]]
    Y.nl=Y[-folds[[l]]]
    
    T.l=T[folds[[l]]]
    T.nl=T[-folds[[l]]]
    
    X.l=X[folds[[l]],]
    X.nl=X[-folds[[l]],]
    
    n.l=length(T.l)
    n.nl=length(T.nl)
    
    # get stage 1 (on nl)
    stage1_estimators<-get_stage1(Y.nl,T.nl,X.nl,
                                  p0,D_LB,D_add,max_iter,dict,
                                  alpha_estimator,gamma_estimator,
                                  c_RR,c_CEF)
    alpha_hat=stage1_estimators[[1]]
    gamma_hat=stage1_estimators[[2]]
    
    # debugging
    print(paste0('fold: ',l))
    
    #get stage 2 (on l) 
    
    #psi_star
    Psi_tilde.l=rep(0,n.l)
    for (i in 1:n.l){
      if(bias){ #plug-in
        Psi_tilde.l[i]=psi_tilde_bias(Y.l[i],T.l[i],X.l[i,],m,alpha_hat,gamma_hat) # without subtracting theta_hat
      }else{ #DML
        Psi_tilde.l[i]=psi_tilde(Y.l[i],T.l[i],X.l[i,],m,alpha_hat,gamma_hat) # without subtracting theta_hat
      }
    }
    
    Psi_tilde=c(Psi_tilde,Psi_tilde.l)
    Y_new=c(Y_new,Y.l)
    T_new=c(T_new,T.l)
    
    #print(paste0(round(mean(Psi_tilde.l),2)))
    
  }
  
  #means
  T_bar=mean(T)
  TY_bar=mean(Y*T)
  
  #point estimation
  est=mean(Psi_tilde)
  ATT=(TY_bar-est)/T_bar
  
  #influences
  Psi=Psi_tilde-est
  
  #V
  v11=mean(Psi^2)
  v12=mean(Psi*T_new*Y_new) #since Psi re-ordered
  v13=mean(Psi*T_new) #since Psi re-ordered
  
  v21=v12
  v22=var(T*Y)
  v23=cov(T*Y,T)
  
  v31=v13
  v32=v23
  v33=var(T)
  
  V=matrix(c(v11,v21,v31,v12,v22,v32,v13,v23,v33),nrow=3)
  
  #H
  h1=-1/T_bar
  h2=1/T_bar
  h3=(est-TY_bar)/(T_bar^2)
  H=matrix(c(h1,h2,h3),nrow=1)
  
  #delta method
  var=H %*% V %*% t(H)
  
  se=sqrt(var/n)
  
  out<-c(table(T)[[2]],table(T)[[1]],ATT,se)
  
  return(out)
}

printer<-function(spec1){
  print(paste(" treated: ",spec1[1], " untreated: ", spec1[2], "   ATT:    ",round(spec1[3],2), "   SE:   ", round(spec1[4],2), sep=""))
}

for_tex<-function(spec1){
  print(paste(" & ",spec1[1], " & ", spec1[2], "   &    ",round(spec1[3],2), "   &   ", round(spec1[4],2), sep=""))
}