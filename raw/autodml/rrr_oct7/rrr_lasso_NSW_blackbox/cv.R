L=5

cv_rr<-function(Y,T,X,p0,D_LB,D_add,max_iter,dict,is_alpha,is_lasso,c){
  n=nrow(X)
  folds <- split(sample(n, n,replace=FALSE), as.factor(1:L))
  
  cv_loss=numeric(0)
  
  for (l in 1:L){
    
    # debugging
    print(paste0('fold: ',l))
    
    Y.l=Y[folds[[l]]]
    Y.nl=Y[-folds[[l]]]
    
    T.l=T[folds[[l]]]
    T.nl=T[-folds[[l]]]
    
    X.l=X[folds[[l]],]
    X.nl=X[-folds[[l]],]
    
    n.l=length(T.l)
    n.nl=length(T.nl)
    
    # get stage 1 (on nl)
    rho.nl=RMD_stable(Y.nl,T.nl,X.nl,p0,D_LB,D_add,max_iter,dict,is_alpha,is_lasso,c)
    
    # get cv loss (on l)
    MNG.l<-get_MNG(Y.l,T.l,X.l,dict)
    M.l=MNG.l[[1]]
    N.l=MNG.l[[2]]
    G.l=MNG.l[[3]]
    if(is_alpha){ # RR
      cv_loss.l=-2 * M.l %*% rho.nl + rho.nl %*% G.l %*% rho.nl
    }else{ # CEF
      cv_loss.l=-2 * N.l %*% rho.nl + rho.nl %*% G.l %*% rho.nl
    }

    # debugging
    #print(paste0('fold: ',l))
    cv_loss=c(cv_loss,cv_loss.l)
    
  }
  
  cv=mean(cv_loss)
  return(cv)
  
}

get_cv_rr<-function(Y,T,X,p0,D_LB,D_add,max_iter,dict,is_alpha,is_lasso,c_vals){
  
  n_vals=length(c_vals)
  cv_vals=rep(NA,n_vals)
  
  for (i in 1:n_vals){
    cv_vals[i]=cv_rr(Y,T,X,p0,D_LB,D_add,max_iter,dict,is_alpha,is_lasso,c_vals[i])
  }
  
  #print(paste0(cv_vals))
  c_star=c_vals[which.min(cv_vals)]
  return(c_star)
}