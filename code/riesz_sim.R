
###########
#Libraries#
###########
rm(list = ls())
library("hdm")
library("foreign")
library("dplyr") 
library("ggplot2")
library("quantreg") #used for rq.fit.sfn
library("nnet")	#used for mulitnom
library("randomForest")
library("keras")
library("data.table")

#################
#Basic Functions#
#################

#Euclidean Norm
two.norm <- function(x){
  return(sqrt(x %*% x))
} 

#dictionary that converts d and z to vector with intercept term
#used to calculate alpha and gamma
b<-function(d,z){
  return(c(1,d,z))
}

#dictionary that converts d and z to vector with intercept term and interaction term
#used to calculate alpha and gamma
b2<-function(d,z){
  return(c(1,d,z,d*z))
}

#identity dictionary
#used with elasticities
b.PD <- function(d,z){
  return(d)
}

#m function mentioned in the introduction for ATE
m.ATE<-function(y,d,z,gamma){ #all data arguments to make interchangeable with m2
  return(gamma(1,z)-gamma(0,z))
}

#m function for finding elasticites
m.PD<-function(y,x,x.up,x.down,delta,gamma){ #all data arguments to make interchangeable with m2
  return((gamma(x.up,0)-gamma(x.down,0))/delta) #0s in gamma are placeholders because gamma requires two arguments
}

#used to calculate matrix N, which is used to calculate gamma_hat and perform regression
m2<-function(y,d,z,gamma){
  return(y*gamma(d,z))
}

#returns number of columns of X
get.p <- function(args,type,X,b){
  if (type == "ATE"){
    return(length(b(args$T[1],X[1,])))
  } else {
    return(ncol(X))
  }
}

#returns number of data points
get.n <- function(args,type,X){
  if (type == "ATE"){
    return(length(args$T))
  } else {
    return(nrow(X))
  }
}

#psi_tilde for debiased machine learning (actual parameter, debiased moment function)
psi_tilde<-function(args,type,y,z,alpha,gamma){
  if (type=="ATE"){
    d = args$T
    #print(d)
    return(m.ATE(y,d,z,gamma)+alpha(d,z)*(y-gamma(d,z)))
  } else if ( type=="ATE_dml" ){
    d = args$T
    r0hat = alpha(d,z)
    alpha_cons = (d/r0hat) - ((1-d)/(1-r0hat))
    return(m.ATE(y,d,z,gamma)+alpha_cons*(y-gamma(d,z)))
  } else {
    x.up = args$X.up
    x.down = args$X.down
    delta = args$delta
    return(m.PD(y,z,x.up,x.down,delta,gamma)+alpha(z,0)*(y-gamma(z,0)))
  }
}
#psi_tilde with a bias. Used to demonstrate shortfall of not accounting for bias. Won't achieve actual debiased machine learning
psi_tilde_bias<-function(args,type,y,z,alpha,gamma){
  if (type=="ATE"){
    d = args$d
    return(m.ATE(y,d,z,gamma))
  } else {
    x.up = args$X.up
    x.down = args$X.down
    delta = args$delta
    return(m.PD(y,z,x.up,x.down,delta,gamma))
  }
}
#Gets M,N and G matrices. These are necessary for calculating the riesz representer, which is discussed in appendix C, and also calculating gamma. Outputs will have dimensions: M (p x 1), N (p x 1), G (p x p)
get_MNG<-function(args, type, Y,X,b){
  
  if (type == "ATE"){
    
    T = args$T
    
    p=length(b(T[1],X[1,]))
    n.nl=length(T)
    
    B=matrix(0,n.nl,p)
    M=matrix(0,p,n.nl)
    N=matrix(0,p,n.nl)
    #Takes average by collecting into a matrix and averaging the rows
    for (i in 1:n.nl){
      B[i,]=b(T[i],X[i,])
      M[,i]=m.ATE(Y[i],T[i],X[i,],b)
      N[,i]=m2(Y[i],T[i],X[i,],b) 
    }
    
    M_hat=rowMeans(M)
    N_hat=rowMeans(N)
    G_hat=t(B)%*%B/n.nl
    
    return(list(M_hat,N_hat,G_hat,B))
  } else {
    
    X.up = args$X.up
    X.down = args$X.down
    delta = args$delta
    #print(X.up)
    #print(X.down)
    p=ncol(X)
    n=nrow(X)
    
    M=matrix(0,p,n)
    N=matrix(0,p,n)
    
    for (i in 1:n){ #simplifications since dictionary b is the identity
      M[,i]=(X.up[i,]-X.down[i,])/delta #since m(w,b)=(x.up-x.down)/delta
      N[,i]=Y[i]*X[i,] #since m2(w,b)=y*x
    }
    
    M_hat=rowMeans(M) #since m(w,b)=dx
    N_hat=rowMeans(N)
    G_hat=t(X)%*%X/n
    
    return(list(M_hat,N_hat,G_hat,X))
  }
}

setwd("/Users/ar8787/Dropbox/data_eco/Udesa/mater_thesis")



get_data<-function(df,spec,quintile,type = "ATE"){
  ###  Inputs
  ###  df: dataset that is being used. In this case a .dta but could be another format
  ###  spec: value of 1,2,or 3 for ATE or 1,2 for PD
  ###  quintile: integer in [0,5]. 0 means all data, 1-5 indicates which quintile of income distribution to filter for
  ### type: what quantity you're trying to estimate (ATE or PD)
  ###  Output
  ###  List of (Y,D,X) or (Y,X,X.up,X.down,delta)
  
  
  if (type == "ATE"){
    df <- na.omit(df)
    Y <- df[,"y"]
    D <- df[,"d"]
    #  All three specifications incorporate the same variables, but differ in the number of higher order terms. X.L has no higher order terms. X.H has higher order terms up to the 8th degree for some variables. X.vH incorporates higher order and interaction terms.
    ## low dim specification
    # X.L <- cbind(df[,"age"], df[,"inc"], df[,"educ"], df[,"fsize"], df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
    
    cols = c('host_listings_count',
             'host_has_profile_pic', 'host_identity_verified', 'Distance_km',
             'accommodation_type_1', 'accommodation_type_2', 'entire_home_apt',
             'private_room', 'accommodates', 'bathrooms', 'bedrooms', 'Real_bed',
             'Wireless_Internet', 'Breakfast', 'Free_Parking', 'instant_bookable',
             'cancellation_policy', 'Smoking_Allowed',
             'require_guest_profile_picture', 'require_guest_phone_verification',
             'Reviews_per_years', 'review_scores_rating')
    X.L <- as.matrix(df[, 3:maxp])
    
    
    # cols = c('host_listings_count',
    # 'host_has_profile_pic', 'host_identity_verified', 'Distance_km',
    # 'accommodation_type_1', 'accommodation_type_2', 'entire_home_apt',
    # 'private_room', 'accommodates', 'bathrooms', 'bedrooms', 'Real_bed',
    # 'Wireless_Internet', 'Breakfast', 'Free_Parking', 'instant_bookable',
    # 'cancellation_policy', 'Smoking_Allowed',
    # 'require_guest_profile_picture', 'require_guest_phone_verification',
    # 'Reviews_per_years', 'review_scores_rating')
    # X.L <- as.matrix(df[, cols])
    # ## high dim specification.
    # cols.h <- c('host_has_profile_pic',
    #             'host_identity_verified',
    #             'accommodation_type_1',
    #             'accommodation_type_2',
    #             'entire_home_apt',
    #             'private_room',
    #             'Real_bed',
    #             'Wireless_Internet',
    #             'Breakfast',
    #             'Free_Parking',
    #             'instant_bookable',
    #             'cancellation_policy',
    #             'Smoking_Allowed',
    #             'require_guest_profile_picture',
    #             'require_guest_phone_verification')
    # X.H <- cbind(poly(df[,"accommodates"], 6, raw=TRUE),
    #              poly(df[,"bathrooms"], 6, raw=TRUE),
    #              poly(df[,"bedrooms"], 6, raw=TRUE),
    #              poly(df[,"Distance_km"], 2, raw=TRUE),
    #              poly(df[,"review_scores_rating"], 4, raw=TRUE),
    #              poly(df[,"Reviews_per_years"], 4, raw=TRUE),
    #              poly(df[,"host_listings_count"], 4, raw=TRUE),
    #              as.matrix(df[, cols.h]))
    # 
    # ## very high dim specification
    # cols.vh <- c('host_has_profile_pic',
    #              'host_identity_verified',
    #              'accommodation_type_1',
    #              'accommodation_type_2',
    #              'entire_home_apt',
    #              'private_room',
    #              'Real_bed',
    #              'Wireless_Internet',
    #              'Breakfast',
    #              'Free_Parking',
    #              'instant_bookable',
    #              'cancellation_policy',
    #              'Smoking_Allowed',
    #              'require_guest_profile_picture',
    #              'require_guest_phone_verification')
    # X.vH=model.matrix(~(poly(df[,"accommodates"], 6, raw=TRUE) +
    #                       poly(df[,"bathrooms"], 6, raw=TRUE) +
    #                       poly(df[,"bedrooms"], 6, raw=TRUE) +
    #                       poly(df[,"Distance_km"], 2, raw=TRUE) +
    #                       poly(df[,"review_scores_rating"], 4, raw=TRUE) +
    #                       poly(df[,"Reviews_per_years"], 4, raw=TRUE) +
    #                       poly(df[,"host_listings_count"], 4, raw=TRUE) +
    #                       df[,"entire_home_apt"] +
    #                       df[,"Breakfast"] +
    #                       df[,"cancellation_policy"] +
    #                       df[,"Smoking_Allowed"] +
    #                       df[,"require_guest_phone_verification"])^2)
    # X.vH=X.vH[,-1]
    
    #Selects for desires specification
    if (spec==1){
      X=X.L
    } else if (spec==2) {
      X=X.H
    } else {
      X=X.vH
    }
    
    X <- scale(X,center=TRUE,scale=TRUE) #Centers and scales X
    n=nrow(X)
    #Impose common support. More discussion in Appendix B.
    p.1 <- multinom(D~X-1, trace=FALSE)$fitted.values
    indexes.to.drop <- which(p.1 < min(p.1[D==1]) | max(p.1[D==1]) < p.1)
    if (length(indexes.to.drop)==0) {indexes.to.drop <- n+1}	#R throws a wobbly if [-indexes.to.drop] is negating an empty set. 
    n.per.treatment <- as.vector(table(D[-indexes.to.drop]))
    n.trim <- n.per.treatment[1]+n.per.treatment[2]
    
    Y.trimmed=Y[-indexes.to.drop]
    D.trimmed=D[-indexes.to.drop]
    X.trimmed=X[-indexes.to.drop,]
    
    if (spec==1){
      inc=X.trimmed[,2]
    } else if (spec==2) {
      inc=X.trimmed[,7]
    } else {
      inc=X.trimmed[,7]
    }
    #Filters for desired quintile
    if (quintile>0){
      q <- ntile(inc, 5)
      Y.q=Y.trimmed[q==quintile]
      D.q=D.trimmed[q==quintile]
      X.q=X.trimmed[q==quintile,]
    } else {
      Y.q=Y.trimmed
      D.q=D.trimmed
      X.q=X.trimmed
    }
    
    return(list(Y.q,D.q,X.q))
    
  } else {
    
    names(df)
    
    N=nrow(df)
    
    # take logs of continuous vars
    cols.to.log<-c("gas","price","income","age","distance")
    df[,cols.to.log]<-lapply(df[,cols.to.log],log)
    
    # output
    Y=df$gas
    
    # construct vars
    df$age2<-((df$age)^2)
    df$income2<-((df$income)^2)
    
    # ensure diff memory location
    df.up<-data.frame(df)
    df.down<-data.frame(df)
    
    # construct df.up and df.down
    prices<-df$price
    delta=sd(prices)/4
    
    df.up$price<-prices+delta/2
    df.down$price<-prices-delta/2
    
    df$price2<-(prices)^2
    df.up$price2<-(df.up$price)^2
    df.down$price2<-(df.down$price)^2
    
    # specification
    if (spec==1){
      formula<- ~ (factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+price+price2+
        price:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year)))+
        price2:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year)))+urban+youngsingle+distance+distance^2+
        age+age2 + income +income2
    } else {
      formula<- ~ (factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+price+price2+
        price:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+age+age2+income+income2)+
        price2:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+age+age2+income+income2)+
        urban+youngsingle+distance+distance^2+
        age+age2 + income +income2
    }
    
    regressors<-model.matrix(formula,data=df)
    regressors.up<-model.matrix(formula, data=df.up)
    regressors.down<-model.matrix(formula, data=df.down)
    
    # check
    dim(regressors)
    dim(regressors.up)
    dim(regressors.down)
    
    if (quintile>0){
      q <- ntile(df$income, 5)
      Y.q=Y[q==quintile]
      regressors.q=regressors[q==quintile,]
      regressors.up.q=regressors.up[q==quintile,]
      regressors.down.q=regressors.down[q==quintile,]
    } else {
      Y.q=Y
      regressors.q=regressors
      regressors.up.q=regressors.up
      regressors.down.q=regressors.down
    }
    
    return(list(Y.q,regressors.q,regressors.up.q,regressors.down.q,delta))  
    
  }
  
}

l=0.1

#Dantzig learning algorithm to tune regularization parameter
RMD_dantzig <- function(M, G, D, lambda=0, sparse = TRUE) {
  
  p <- ncol(G)
  zp <- rep(0, p)
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,T,X,TX)
  
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
#Lasso learning algorithm to tune the regularization paremeter of our model
RMD_lasso <- function(M, G, D, lambda=0, 
                      control = list(maxIter = 1000, optTol = 10^(-5), 
                                     zeroThreshold = 10^(-6)), beta.start = NULL) {
  
  p <- ncol(G)
  
  Gt<-G
  Mt<-M
  
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,...)
  lambda_vec=lambda*L*D #v3: insert D here
  
  if (is.null(beta.start)) {
    beta <- rep(0,p) #vs low-dimensional initialization
  }
  else {
    beta <- beta.start
  }
  wp <- beta
  mm <- 1
  while (mm < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      rho=Mt[j]-Gt[j,]%*%beta+Gt[j,j]*beta[j]
      z=Gt[j,j]
      
      if (sum(is.na(rho)) >= 1) {
        beta[j] <- 0
        next
      }
      if (rho < -1 * lambda_vec[j]) 
        beta[j] <- (rho+lambda_vec[j])/z
      if (abs(rho) <= lambda_vec[j]) 
        beta[j] <- 0
      if (rho > lambda_vec[j]) 
        beta[j] <- (rho-lambda_vec[j])/z
    }
    wp <- cbind(wp, beta)
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    mm <- mm + 1
  }
  w <- beta
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = mm))
}

#Performs theoretical iteration to help tune our regularization parameter
get_D <- function(args,type,Y,X,m,rho_hat,b){
  #n=nrow(X)
  #p=length(b(T[1],X[1,]))
  n = get.n(args,type,X)
  p = get.p(args,type,X,b)
  
  df=matrix(0,p,n)
  if (type== "ATE"){
    T = args$T
    for (i in 1:n){
      df[,i]=b(T[i],X[i,])*as.vector(rho_hat %*% b(T[i],X[i,]))-m.ATE(Y[i],T[i],X[i,],b)
    }
  } else {
    X.up = args$X.up
    X.down = args$X.down
    delta = args$delta
    for (i in 1:n){
      df[,i]=X[i,]*as.vector(rho_hat %*% X[i,])-m.PD(Y[i],X[i,],X.up[i,],X.down[i,],delta,b)
    }
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

c=0.5
alpha=0.1
tol=1e-6
#Calculates ML function for alpha and gamma
RMD_stable<-function(args,type,Y,X,p0,D_LB,D_add,max_iter,b,is_alpha = TRUE,is_lasso = TRUE){
  ###Inputs
  ###args: list of differing arguments between ATE and PD
  ###   ATE -> args = list(T)
  ###   PD -> args = list(X.up,X.down,delta)
  ###Y = outcome variable
  ###X = input variables
  ###type = ATE or PD
  ###p0: initial parameter for p
  ###D_LB: Lower bound on D
  ###D_add: value added to D
  ###max_iter: maximum iterations of chosen ML algorithm
  ###b: dictionary 
  ###is_alpha: True = return alpha, False = return gamma
  ###is_lasso: True = use lasso, False = use dantzig
  ###Output
  ###rho_hat or beta_hat, a vector to estimate alpha or gamma rescpectively by taking the dot product with b
  k=1
  
  #p=length(b(T[1],X[1,]))
  #n=length(T)
  
  p = get.p(args,type,X,b)
  n = get.n(args,type,X)
  
  # low-dimensional moments
  X0=X[,1:p0]
  if (type == "ATE"){
    args0 = list(T = args$T)
  } else {
    args0 = list(X.up = args$X.up[,1:p0], X.down = args$X.down[,1:p0], delta = args$delta)
  }
  MNG0<-get_MNG(args0, type, Y,X0,b)
  M_hat0=MNG0[[1]]
  N_hat0=MNG0[[2]]
  G_hat0=MNG0[[3]]
  
  # initial estimate
  rho_hat0=solve(G_hat0,M_hat0)
  rho_hat=c(rho_hat0,rep(0,p-ncol(G_hat0)))
  beta_hat0=solve(G_hat0,N_hat0)
  beta_hat=c(beta_hat0,rep(0,p-ncol(G_hat0)))
  
  # moments
  MNG<-get_MNG(args,type,Y,X,b)
  M_hat=MNG[[1]]
  N_hat=MNG[[2]]
  G_hat=MNG[[3]]
  
  # penalty
  lambda=c*qnorm(1-alpha/(2*p))/sqrt(n) # snippet
  
  if(is_alpha){ 
    ###########
    # alpha_hat
    ###########
    diff_rho=1
    #Loop through max_iter times or until the change in rho between iterations is less than tol
    while(diff_rho>tol & k<=max_iter){
      
      # previous values
      rho_hat_old=rho_hat+0
      
      # normalization
      D_hat_rho=get_D(args,type,Y,X,m,rho_hat_old,b)
      D_hat_rho=pmax(D_LB,D_hat_rho)
      D_hat_rho=D_hat_rho+D_add
      
      # RMD estimate
      if(is_lasso){
        rho_hat=RMD_lasso(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }else{
        rho_hat=RMD_dantzig(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }
      
      # difference
      diff_rho=two.norm(rho_hat-rho_hat_old)
      k=k+1
      
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(rho_hat)
    
  } else { 
    ###########
    # gamma_hat
    ###########
    diff_beta=1
    #Loop through max_iter times or until the change in rho between iterations is less than tol
    while(diff_beta>tol & k<=max_iter){
      
      # previous values
      beta_hat_old=beta_hat+0
      
      # normalization
      D_hat_beta=get_D(args,type,Y,X,m2,beta_hat_old,b)
      D_hat_beta=pmax(D_LB,D_hat_beta)
      D_hat_beta=D_hat_beta+D_add
      
      # RMD estimate
      if(is_lasso){
        beta_hat=RMD_lasso(N_hat, G_hat, D_hat_beta, lambda)$coefficients
      }else{
        beta_hat=RMD_dantzig(N_hat, G_hat, D_hat_beta, lambda)$coefficients
      }
      
      # difference
      diff_beta=two.norm(beta_hat-beta_hat_old)
      k=k+1
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(beta_hat)
    
  }
}

arg_Forest<- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
arg_Nnet<- list(size=8,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)


get_stage1<-function(args,type,Y,X,p0,D_LB,D_add,max_iter,b,alpha_estimator = 1,gamma_estimator = 1){
  ###Inputs
  ###args: list of differing arguments between ATE and PD
  ###   ATE -> args = list(T)
  ###   PD -> args = list(X.up,X.down,delta)
  ###Y = outcome variable
  ###X = input variables
  ###type = ATE or PD
  ###p0: initial parameter for p
  ###D_LB: Lower bound on D
  ###D_add: value added to D
  ###max_iter: maximum iterations of chosen ML algorithm
  ###b: dictionary 
  ###alpha_estimator: numerical value indicating which model to use for the alpha estimator
  ###0 dantzig, 1 lasso
  ###gamma_estimator: numerical value indicating which model to use for the gamma estimator
  ###0 dantzig, 1 lasso, 2 rf, 3 nn
  ###Output
  ###alpha_hat and gamma_hat estimators
  
  #For ATE
  #p=length(b(T[1],X[1,]))
  #n=length(T)
  
  #For PD
  #p = dim(X)[2]
  #n = dim(X)[1]
  
  p = get.p(args,type,X,b)
  n = get.n(args,type,X)
  MNG<-get_MNG(args,type,Y,X,b)
  B=MNG[[4]]
  
  ###########
  # alpha hat
  ###########
  if(alpha_estimator==0){ # dantzig
    
    rho_hat=RMD_stable(args,type,Y,X,p0,D_LB,D_add,max_iter,b,1,0)
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  } else if(alpha_estimator==1){ # lasso
    
    rho_hat=RMD_stable(args,type,Y,X,p0,D_LB,D_add,max_iter,b,1,1)
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  }
  
  ###########
  # gamma hat
  ###########
  if(gamma_estimator==0){ # dantzig
    
    beta_hat=RMD_stable(args,type,Y,X,p0,D_LB,D_add,max_iter,b,0,0)
    gamma_hat<-function(d,z){
      return(b(d,z)%*%beta_hat)
    }
    
  } else if(gamma_estimator==1){ # lasso
    
    beta_hat=RMD_stable(args,type,Y,X,p0,D_LB,D_add,max_iter,b,0,1)
    gamma_hat<-function(d,z){ 
      return(b(d,z)%*%beta_hat)
    }
    
  } else if(gamma_estimator==2){ # random forest
    
    forest<- do.call(randomForest, append(list(x=B,y=Y), arg_Forest))
    gamma_hat<-function(d,z){
      return(predict(forest,newdata=b(d,z), type="response"))
    }
    
  } else if(gamma_estimator==3){ # neural net
    
    # scale down, de-mean, run NN, scale up, remean so that NN works well
    maxs_B <- apply(B, 2, max)
    mins_B <- apply(B, 2, min)
    
    maxs_Y<-max(Y)
    mins_Y<-min(Y)
    
    # hack to ensure that constant covariates do not become NA in the scaling
    const=maxs_B==mins_B
    keep=(1-const)*1:length(const)
    
    NN_B<-B
    NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maxs_B[keep] - mins_B[keep])
    
    NN_Y<-scale(Y, center = mins_Y, scale = maxs_Y - mins_Y)
    
    nn<- do.call(nnet, append(list(x=NN_B,y=NN_Y), arg_Nnet))
    gamma_hat<-function(d,z){
      
      test<-t(as.vector(b2(d,z)))
      NN_b<-test
      NN_b[,keep]<-scale(t(NN_b[,keep]), 
                         center = mins_B[keep], 
                         scale = maxs_B[keep] - mins_B[keep])
      
      NN_Y_hat<-predict(nn,newdata=NN_b)
      Y_hat=NN_Y_hat*(maxs_Y-mins_Y)+mins_Y
      
      return(Y_hat)
    }
    
  } else if(gamma_estimator==4){ # 2 layer NN (keras)
    
    
    # scale down, de-mean, run NN, scale up, remean so that NN works well
    # hack to ensure that constant covariates do not become NA in the scaling
    maxs_B <- apply(B, 2, max)
    mins_B <- apply(B, 2, min)
    maxs_Y<-max(Y)
    mins_Y<-min(Y)
    const=maxs_B==mins_B
    keep=(1-const)*1:length(const)
    NN_B<-B
    NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maxs_B[keep] - mins_B[keep])
    NN_Y<-scale(Y, center = mins_Y, scale = maxs_Y - mins_Y)
    
    # Corrected build_model function
    build_model <- function() {
      # Input shape should be the number of features (columns in NN_B)
      input <- layer_input(shape = ncol(NN_B))
      
      # Define the model architecture using pipe operator
      output <- input %>%
        layer_dense(units = 8, activation = "relu") %>%
        layer_dense(units = 8, activation = "relu") %>%
        layer_dense(units = 1)
      
      # Create and compile the model
      model <- keras_model(inputs = input, outputs = output)
      
      model %>% compile(
        optimizer = "rmsprop",
        loss = "mse",
        metrics = "mae"
      )
      
      return(model)
    }
    
    
    
    
    
    
    # use package
    model <- build_model()
    num_epochs <- 100
    model %>% fit(NN_B, NN_Y,
                  epochs = num_epochs, batch_size = 1, verbose = 0)
    
    gamma_hat<-function(d,z){
      
      NN_b<-t(as.vector(b(d,z))) # test point
      NN_b[,keep]<-scale(t(NN_b[,keep]), 
                         center = mins_B[keep], 
                         scale = maxs_B[keep] - mins_B[keep])
      
      # 2 layer NN (keras)
      NN_Y_hat<-model %>% predict(NN_b, verbose = 0)
      
      Y_hat=NN_Y_hat*(maxs_Y-mins_Y)+mins_Y
      
      return(Y_hat)
    }
    
  }
  
  # X_const <- cbind(Intercept = 1, X)
  # t.coef <- cv.glmnet(X_const, args$T, family="binomial", alpha=1)
  # prob.t<-function(d,z){ 
  #   # I am adding the intercept
  #   return(predict(t.coef,c(1,z)))
  # } 
  
  t.coef <- rlasso(X, args$T,
                   intercept = TRUE,
                   penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, lambda.start = NULL, c = 1.1))$coefficients
  prob.t<-function(d,z){
    # I am adding the intercept
    return( c(1,z) %*% t.coef)
  }
  
  li_hat=RMD_stable(args,type,args$T,X,p0,D_LB,D_add,max_iter,b,0,1)
  alpha.aproxx<-function(d,z){ 
    return(b(d,z)%*%li_hat)
  }
  
  
  return(list(alpha_hat,gamma_hat, prob.t, alpha.aproxx))
  
}

L=5 #parameter for number of subsets to split data into

rrr<-function(args,type,Y,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias){
  if (type=='ATE'){
    T = args$T
  } else {
    X.up = args$X.up
    X.down = args$X.down
    delta = args$delta
  }
  n=nrow(X)
  #Split up into I_1,...I_n
  folds <- split(sample(n, n,replace=FALSE), as.factor(1:L))
  
  Psi_tilde.new=numeric(0)
  Psi_tilde.dml=numeric(0)
  Psi_tilde.adml=numeric(0) #Initialize Psi_tilde as 0. This will become a running list of
  #Iterating through I_1,,,,I_n
  for (l in 1:L){
    #Breaks into two sets. variables with .l are not in X_i used for training. variables with .nl are used for calculating theta
    Y.l=Y[folds[[l]]]
    Y.nl=Y[-folds[[l]]]
    
    X.l=X[folds[[l]],]
    X.nl=X[-folds[[l]],]
    
    
    
    if (type=="ATE"){
      
      T.l=T[folds[[l]]]
      T.nl=T[-folds[[l]]]
      
      n.l=length(T.l)
      n.nl=length(T.nl)
      
      stage1.args = list(T = T.nl)
      
      
    } else {
      X.up.l=X.up[folds[[l]],]
      X.up.nl=X.up[-folds[[l]],]
      
      X.down.l=X.down[folds[[l]],]
      X.down.nl=X.down[-folds[[l]],]
      
      n.l=nrow(X.l)
      n.nl=nrow(X.nl)
      
      stage1.args = list(X.up = X.up.nl, X.down = X.down.nl, delta = delta)
    } 
    
    # get stage 1 (on nl)
    stage1_estimators<-get_stage1(stage1.args, type, Y.nl,X.nl,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator)
    alpha_hat=stage1_estimators[[1]]
    gamma_hat=stage1_estimators[[2]]
    prob.t=stage1_estimators[[3]]
    alpha.aproxx=stage1_estimators[[4]]
    print(paste0('fold: ',l))
    
    #get stage 2 (on l)
    #psi_star
    Psi_tilde.l.adml=rep(0,n.l)
    for (i in 1:n.l){
      if (type=="ATE"){psi.args = list(T = T.l[i])} 
      else {psi.args = list(X.up = X.up.l[i,],X.down = X.down.l[i,],delta = delta)}
      if(bias){ #plug-in
        Psi_tilde.l.adml[i]=psi_tilde_bias(psi.args, type, Y.l[i],X.l[i,],alpha_hat,gamma_hat) # without subtracting theta_hat
      }else{ #DML
        Psi_tilde.l.adml[i]=psi_tilde(psi.args, type, Y.l[i],X.l[i,],alpha_hat,gamma_hat) # without subtracting theta_hat
      }
    }
    
    
    Psi_tilde.l.dml=rep(0,n.l)
    for (i in 1:n.l){
      if (type=="ATE"){psi.args = list(T = T.l[i])} 
      else {psi.args = list(X.up = X.up.l[i,],X.down = X.down.l[i,],delta = delta)}
      if(bias){ #plug-in
        Psi_tilde.l.dml[i]=psi_tilde_bias(psi.args, "ATE_dml", Y.l[i],X.l[i,],prob.t,gamma_hat) # without subtracting theta_hat
      }else{ #DML
        Psi_tilde.l.dml[i]=psi_tilde(psi.args, "ATE_dml", Y.l[i],X.l[i,],prob.t,gamma_hat) # without subtracting theta_hat
      }
    }
    
    
    Psi_tilde.l.new=rep(0,n.l)
    for (i in 1:n.l){
      if (type=="ATE"){psi.args = list(T = T.l[i])} 
      else {psi.args = list(X.up = X.up.l[i,],X.down = X.down.l[i,],delta = delta)}
      if(bias){ #plug-in
        Psi_tilde.l.new[i]=psi_tilde_bias(psi.args, type, Y.l[i],X.l[i,],alpha.aproxx,gamma_hat) # without subtracting theta_hat
      }else{ #DML
        Psi_tilde.l.new[i]=psi_tilde(psi.args, type, Y.l[i],X.l[i,],alpha.aproxx,gamma_hat) # without subtracting theta_hat
      }
    }
    
    
    Psi_tilde.adml=c(Psi_tilde.adml,Psi_tilde.l.adml)
    Psi_tilde.dml=c(Psi_tilde.dml,Psi_tilde.l.dml)
    Psi_tilde.new=c(Psi_tilde.new,Psi_tilde.l.new)
    
    
    
  }
  
  #point estimation
  ate.adml=mean(Psi_tilde.adml)
  #influences
  Psi.adml=Psi_tilde.adml-ate.adml
  var.adml=mean(Psi.adml^2)
  se.adml=sqrt(var.adml/n)
  
  ate.dml=mean(Psi_tilde.dml, na.rm = TRUE)
  Psi.dml=Psi_tilde.dml-ate.dml
  var.dml=mean(Psi.dml^2, na.rm = TRUE)
  se.dml=sqrt(var.dml/n)
  
  
  ate.new=mean(Psi_tilde.new, na.rm = TRUE)
  Psi.new=Psi_tilde.new-ate.new
  var.new=mean(Psi.new^2, na.rm = TRUE)
  se.new=sqrt(var.new/n)
  
  #Returns ATE and SE
  if (type == "ATE"){
    out<-c(table(T)[[2]],table(T)[[1]], 
           'ate.adml'= ate.adml, 'se.adml' = se.adml, 
           'ate.dml' = ate.dml, 'se.dml' = se.dml, 
           'ate.new' = ate.new, 'se.dml' = se.new)
  } else {
    out<-c(n,ate,se)
  }
  
  return(out)
}
#Prints out the results of rrr in a way that is easy to read
printer<-function(spec1,type){
  if (type == "ATE"){
    print(paste(" treated: ",spec1[1], " untreated: ", spec1[2], "   ATE:    ",round(spec1[3],2), "   SE:   ", round(spec1[4],2), sep=""))
  } else {
    print(paste("   n:    ",spec1[1],"   AD:    ",round(spec1[2],2), "   SE:   ", round(spec1[3],2), sep=""))
  }
}
#Prints results of rrr in format suitable for LaTex
for_tex<-function(spec1,type){
  if (type == "ATE"){
    print(paste(" & ",spec1[1], " & ", spec1[2], "   &    ",round(spec1[3],2), "   &   ", round(spec1[4],2), sep=""))
  } else {
    print(paste(" & ",spec1[1]," & ",round(spec1[2],2), "   &   ", round(spec1[3],2), sep=""))
  }
}

type = "ATE" #either ATE for average treatment effect or PD for elasticity
max_iter = 1000
for (quintile in 0:0){ #change "0:5" to "0:0" if you just want to run once on whole data set
  
  quintile=0 #comment out this line to iterate through quintiles
  
  print(paste0('quintile: '))
  print(paste0(quintile))
  
  #This if statement sets up data for either ATE or PD
  #note: each clause has its own 'spec' and 'dict' variables, so
  #if you want to change these, make sure you're changing it in 
  #the correct clause
  #Conversely, all other arguments (such as alpha) are outside these clauses,
  #so changing them for one type will change them for either type
  if (type=="ATE"){ 
    df  <- fread("/Users/ar8787/Dropbox/data_eco/Udesa/mater_thesis/data/simdata.csv") %>% 
      as.data.frame() %>% 
      na.omit()
    # df <- df[, 1:11]
    maxp = ncol(df)
    #df = read.csv("Voting.csv")
    
    spec=1 #spec in (1-3)
    quintile=0 #quintile in (1-5). 0 means all quintiles
    data<-get_data(df,spec,quintile) #trimming like Farrell; different than Chernozhukov et al. 
    
    Y=data[[1]]
    T=data[[2]]
    X=data[[3]] #no intercept
    
    
    ##################
    # helper functions
    ##################
    
    
    # dictionary for calculating alpha (discussion in appendix C)
    dict=b # b for partially linear model, b2 for interacted model. note that b2 appears in stage1.R for NN
    p=length(b(T[1],X[1,]))
    args = list(T = T)
  } else {
    df  <- read.dta("gasoline_final_tf1.dta")
    spec = 1 #either 1 or 2
    data<-get_data(df,spec,quintile,type = "PD") #like Chernozhukov and Semenova
    Y=data[[1]]
    X=data[[2]]
    X.up=data[[3]]
    X.down=data[[4]]
    delta=data[[5]]
    
    #test<-get_MNG(Y,X,X.up,X.down,delta)
    #M_hat=test[[1]]
    #N_hat=test[[2]]
    #as.numeric(M_hat)
    #as.numeric(N_hat)
    
    n=nrow(X)
    p=ncol(X)
    dict = b.PD
    args = list(X.up = X.up, X.down = X.down, delta = delta)
    
  }
  
  
  #p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
  p0=ceiling(p/4) 
  if (p>60){
    p0=ceiling(p/40)
    
  }
  
  
  D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
  D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
  max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights
  
  ###########
  # algorithm
  ###########
  
  set.seed(1) # for sample splitting
  alpha_estimator=1
  gamma_estimator=1
  bias=0 #0 for unbiased or 1 for biased
  #alpha_estimator: 0 dantzig, 1 lasso
  #gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn, 4 2-layer nn (2-layer nn only works for ATE)
  results<-rrr(args,type,Y,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias)
  printer(results,type)
  
  
}


set.seed(1) # Random Fores
alpha_estimator=1
gamma_estimator=2
bias=0
results2<-rrr(args,type,Y,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias)
# -0.0348 - (0.0052)


set.seed(1) # for sample splitting
alpha_estimator=1
gamma_estimator=3
bias=0
results3<-rrr(args,type,Y,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias)
# -0.0357 - (0.0056)

# very time consuming
set.seed(1) # for sample splitting
alpha_estimator=1
gamma_estimator=4
bias=0
results4<-rrr(args,type,Y,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias)
# 0.0275407403116969 (0.0059)




folds_list <- results[[2]]
results[[2]]

# Convert each fold (assumed to be a vector) into a character string of indices
folds_df <- data.frame(
  Fold = seq_along(folds_list),
  Indices = sapply(folds_list, function(x) paste(x-1, collapse = ","))
)


# 
# df  <- read.dta("sipp1991.dta")
# Y <- df[,"net_tfa"]
# D <- df[,"e401"]
# X.L <- cbind(df[,"age"], df[,"inc"], df[,"educ"], df[,"fsize"], df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
# 
# X <- scale(X.L,center=TRUE,scale=TRUE) #Centers and scales X
# #Impose common support. More discussion in Appendix B.
# p.1 <- multinom(D~X-1, trace=FALSE)$fitted.values
# indexes.to.drop <- which(p.1 < min(p.1[D==1]) | max(p.1[D==1]) < p.1)
# 
# probs = data.frame("treatment" = append(rep("Treated",length(p.1[D==1])),rep("Untreated",length(p.1[D==0]))),"value" = append(p.1[D==1],p.1[D==0]))  
# 
# ggplot(probs,aes(x = value, fill = treatment))+
#   geom_histogram(alpha = 0.2, col = "black")+
#   geom_vline(xintercept = min(p.1[D==1]),linetype = "dotted")+
#   geom_vline(xintercept = max(p.1[D==1]),linetype = "dotted")
# 
# probs %>%
#   filter(value<=0.25)%>%
#   ggplot(aes(x = value, fill = treatment))+
#   geom_histogram(alpha = 0.2, col = "black")+
#   geom_vline(xintercept = min(p.1[D==1]))
# 
