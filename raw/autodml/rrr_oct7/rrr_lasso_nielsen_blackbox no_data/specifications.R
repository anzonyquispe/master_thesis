get_data<-function(df,milk_indicator,poly_degree,centering){
  
  if(milk_indicator){
    Y<-df$milk_s
    
    # notation of Woodlridge: X for direct effects, Z for correlated random effects
    X=model.matrix(~(poly(milk,poly_degree,raw=T)+poly(expend,poly_degree, raw=T))^4+
                     bread+butter+cereal+chips+coffee+cookies+eggs+ic+oj+salad+soda+soup+water+yogurt+
                     hh+t+n,data=df) #drops NA rows
    
    # dX matrix
    dX<-matrix(0,dim(X)[1],dim(X)[2])
    colnames(dX)<-colnames(X)
    
    # find indices
    milk1.idx<-grep("poly(milk, poly_degree, raw = T)1",colnames(X),fixed=TRUE)
    milk2.idx<-grep("poly(milk, poly_degree, raw = T)2",colnames(X),fixed=TRUE)
    milk3.idx<-grep("poly(milk, poly_degree, raw = T)3",colnames(X),fixed=TRUE)
    milk4.idx<-grep("poly(milk, poly_degree, raw = T)4",colnames(X),fixed=TRUE)
    
    # take derivatives
    dX[,milk1.idx]<-X[,milk1.idx]/df$milk
    dX[,milk2.idx]<-2*X[,milk2.idx]/df$milk
    dX[,milk3.idx]<-3*X[,milk3.idx]/df$milk
    dX[,milk4.idx]<-4*X[,milk4.idx]/df$milk
    
  }else{
    Y<-df$soda_s
    
    # notation of Woodlridge: X for direct effects, Z for correlated random effects
    X=model.matrix(~(poly(soda,poly_degree,raw=T)+poly(expend,poly_degree, raw=T))^4+
                     bread+butter+cereal+chips+coffee+cookies+eggs+ic+milk+oj+salad+soup+water+yogurt+
                     hh+t+n,data=df) #drops NA rows
    
    # dX matrix
    dX<-matrix(0,dim(X)[1],dim(X)[2])
    colnames(dX)<-colnames(X)
    
    # find indices
    soda1.idx<-grep("poly(soda, poly_degree, raw = T)1",colnames(X),fixed=TRUE)
    soda2.idx<-grep("poly(soda, poly_degree, raw = T)2",colnames(X),fixed=TRUE)
    soda3.idx<-grep("poly(soda, poly_degree, raw = T)3",colnames(X),fixed=TRUE)
    soda4.idx<-grep("poly(soda, poly_degree, raw = T)4",colnames(X),fixed=TRUE)
    
    # take derivatives
    dX[,soda1.idx]<-X[,soda1.idx]/df$soda
    dX[,soda2.idx]<-2*X[,soda2.idx]/df$soda
    dX[,soda3.idx]<-3*X[,soda3.idx]/df$soda
    dX[,soda4.idx]<-4*X[,soda4.idx]/df$soda
    
  }

  # preserve hh, t, n
  hh.idx<-grep("hh",colnames(X))
  t.idx<-grep("^t$",colnames(X))
  n.idx<-grep("^n$",colnames(X))
  dX[,hh.idx]<-X[,hh.idx]
  dX[,t.idx]<-X[,t.idx]
  dX[,n.idx]<-X[,n.idx]
  
  # summaries
  X_bar<-as.tibble(X) %>% group_by(hh) %>% summarise_all(.funs = c(mean="mean"))
  X_hat<-X_bar %>% summarise_all(.funs = c(mean="mean"))
  dX_bar<-as.tibble(dX) %>% group_by(hh) %>% summarise_all(.funs = c(mean="mean")) #currently unused. may appear later?
  dX_hat<-dX_bar %>% summarise_all(.funs = c(mean="mean"))
  
  # regressors
  X<-as.tibble(X) %>% dplyr::select(-c(hh,t,n))
  X<-as.matrix(X[,-1])
  Z<-X_bar %>% uncount(n_mean) %>% dplyr::select(-c(hh,t_mean))
  Z<-as.matrix(Z[,-1])
  mu<-X_hat %>% slice(rep(row_number(), length(df$hh))) %>% dplyr::select(-c(hh_mean,t_mean_mean,n_mean_mean))
  mu<-as.matrix(mu[,-1])
  
  # center and scale indirect components
  regressors3<-Z-mu
  regressors4<-kr(Z-mu,X)
  if(centering){
    regressors3<-scale(regressors3)
    regressors4<-scale(regressors4)
  }
  
  regressors<-cbind(1, X,regressors3,regressors4) #intercept
  
  # dregressors; only dX nonzero
  dX<-as.tibble(dX) %>% dplyr::select(-c(hh,t,n))
  dX<-as.matrix(dX[,-1])
  dZ<-matrix(0,dim(dX)[1],dim(dX)[2])
  dmu<-matrix(0,dim(dX)[1],dim(dX)[2])
  dregressors<-cbind(0, dX, dZ-dmu,kr(dZ-dmu,dX)) #intercept
  
  return(list(df$hh,df$t,df$n,Y,regressors,dregressors))
  
}


