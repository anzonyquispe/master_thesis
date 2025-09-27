rm(list = ls())
setwd("/Users/ar8787/Dropbox/data_eco/Udesa/mater_thesis/raw/18515_Data_and_Programs/tutorial")
# Y is net_tfa
# D is e401 - Treatment
# X is covs (low dimmension spec == 1, 2 high, 3, very high)
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
RMD_lasso <- function(M, G, D, lambda=0, control = list(maxIter = 1000, optTol = 10^(-5), 
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


#Euclidean Norm
two.norm <- function(x){
  return(sqrt(x %*% x))
} 

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
  } else {
    x.up = args$X.up
    x.down = args$X.down
    delta = args$delta
    return(m.PD(y,z,x.up,x.down,delta,gamma)+alpha(z,0)*(y-gamma(z,0)))
  }
}
get_data<-function(df,spec,quintile,type = "ATE"){
  ###  Inputs
  ###  df: dataset that is being used. In this case a .dta but could be another format
  ###  spec: value of 1,2,or 3 for ATE or 1,2 for PD
  ###  quintile: integer in [0,5]. 0 means all data, 1-5 indicates which quintile of income distribution to filter for
  ### type: what quantity you're trying to estimate (ATE or PD)
  ###  Output
  ###  List of (Y,D,X) or (Y,X,X.up,X.down,delta)
  
  
  if (type == "ATE"){
    Y <- df[,"net_tfa"]
    D <- df[,"e401"]
    #  All three specifications incorporate the same variables, but differ in the number of higher order terms. X.L has no higher order terms. X.H has higher order terms up to the 8th degree for some variables. X.vH incorporates higher order and interaction terms.
    ## low dim specification
    X.L <- cbind(df[,"age"], df[,"inc"], df[,"educ"], df[,"fsize"], df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
    
    ## high dim specification.
    X.H <- cbind(poly(df[,"age"], 6, raw=TRUE),
                 poly(df[,"inc"], 8, raw=TRUE),
                 poly(df[,"educ"], 4, raw=TRUE),
                 poly(df[,"fsize"], 2, raw=TRUE),
                 df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
    
    ## very high dim specification
    X.vH=model.matrix(~(poly(df[,"age"], 6, raw=TRUE) +
                          poly(df[,"inc"], 8, raw=TRUE) +
                          poly(df[,"educ"], 4, raw=TRUE) +
                          poly(df[,"fsize"], 2, raw=TRUE) +
                          df[,"marr"] +
                          df[,"twoearn"] +
                          df[,"db"] +
                          df[,"pira"] +
                          df[,"hown"])^2)
    X.vH=X.vH[,-1]
    
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

df  <- read.dta("sipp1991.dta")
#df = read.csv("Voting.csv")

spec=1 #spec in (1-3)
#quintile=0 #quintile in (1-5). 0 means all quintiles
quintile <- 0
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
# p number of covariables (9) including the intercept (1) and treatment (1) == 11
args = list(T = T)


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
c=0.5
alpha=0.1
tol=1e-6
set.seed(1) # for sample splitting
type = "ATE"
alpha_estimator=1
gamma_estimator=1
bias=0 #0 for unbiased or 1 for biased
L = 5 # Number of Folds
#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn, 4 2-layer nn (2-layer nn only works for ATE)
# results<-rrr(args,type,Y,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias)
get_D <- function(args,type,Y,X,m,rho_hat,b){
  #n=nrow(X)
  #p=length(b(T[1],X[1,]))
  n = get.n(args,type,X)
  p = get.p(args,type,X,b)
  
  df=matrix(0,p,n)
  if (type== "ATE"){
    T = args$T
    for (i in 1:n){
      df[,i]=b(T[i],X[i,])*as.vector(rho_hat %*% b(T[i],X[i,])) -m.ATE(Y[i],T[i],X[i,],b)
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
  folds <- split(sample(n, n,replace=FALSE), as.factor(1:L)) # sampling with replacement all the folds with the same n of observations
  
  Psi_tilde=numeric(0) #Initialize Psi_tilde as 0. This will become a running list of
  l = 1
  #Iterating through I_1,,,,I_n
  for (l in 1:L){
    #Breaks into two sets. variables with .l are not in X_i used for training. variables with .nl are used for calculating theta
    Y.l=Y[folds[[l]]] # data to generate the parameter
    Y.nl=Y[-folds[[l]]] # data use to train the model alpha and gamma
    
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
    
    print(paste0('fold: ',l))
    
    #get stage 2 (on l)
    #psi_star
    Psi_tilde.l=rep(0,n.l)
    for (i in 1:n.l){
      if (type=="ATE"){psi.args = list(T = T.l[i])} 
      else {psi.args = list(X.up = X.up.l[i,],X.down = X.down.l[i,],delta = delta)}
      if(bias){ #plug-in
        Psi_tilde.l[i]=psi_tilde_bias(psi.args, type, Y.l[i],X.l[i,],alpha_hat,gamma_hat) # without subtracting theta_hat
      }else{ #DML
        Psi_tilde.l[i]=psi_tilde(psi.args, type, Y.l[i],X.l[i,],alpha_hat,gamma_hat) # without subtracting theta_hat
      }
    }
    
    Psi_tilde=c(Psi_tilde,Psi_tilde.l)
    
    
  }
  
  #point estimation
  ate=mean(Psi_tilde)
  
  #influences
  Psi=Psi_tilde-ate
  
  var=mean(Psi^2)
  se=sqrt(var/n)
  #Returns ATE and SE
  if (type == "ATE"){
    out<-c(table(T)[[2]],table(T)[[1]],ate,se)
  } else {
    out<-c(n,ate,se)
  }
  
  return(out)
}









rm(list = ls())
setwd("/Users/ar8787/Dropbox/data_eco/Udesa/mater_thesis/raw/18515_Data_and_Programs/tutorial")
# Y is net_tfa
# D is e401 - Treatment
# X is covs (low dimmension spec == 1, 2 high, 3, very high)

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
  } else {
    x.up = args$X.up
    x.down = args$X.down
    delta = args$delta
    return(m.PD(y,z,x.up,x.down,delta,gamma)+alpha(z,0)*(y-gamma(z,0)))
  }
}
get_data<-function(df,spec,quintile,type = "ATE"){
  ###  Inputs
  ###  df: dataset that is being used. In this case a .dta but could be another format
  ###  spec: value of 1,2,or 3 for ATE or 1,2 for PD
  ###  quintile: integer in [0,5]. 0 means all data, 1-5 indicates which quintile of income distribution to filter for
  ### type: what quantity you're trying to estimate (ATE or PD)
  ###  Output
  ###  List of (Y,D,X) or (Y,X,X.up,X.down,delta)
  
  
  if (type == "ATE"){
    Y <- df[,"net_tfa"]
    D <- df[,"e401"]
    #  All three specifications incorporate the same variables, but differ in the number of higher order terms. X.L has no higher order terms. X.H has higher order terms up to the 8th degree for some variables. X.vH incorporates higher order and interaction terms.
    ## low dim specification
    X.L <- cbind(df[,"age"], df[,"inc"], df[,"educ"], df[,"fsize"], df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
    
    ## high dim specification.
    X.H <- cbind(poly(df[,"age"], 6, raw=TRUE),
                 poly(df[,"inc"], 8, raw=TRUE),
                 poly(df[,"educ"], 4, raw=TRUE),
                 poly(df[,"fsize"], 2, raw=TRUE),
                 df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
    
    ## very high dim specification
    X.vH=model.matrix(~(poly(df[,"age"], 6, raw=TRUE) +
                          poly(df[,"inc"], 8, raw=TRUE) +
                          poly(df[,"educ"], 4, raw=TRUE) +
                          poly(df[,"fsize"], 2, raw=TRUE) +
                          df[,"marr"] +
                          df[,"twoearn"] +
                          df[,"db"] +
                          df[,"pira"] +
                          df[,"hown"])^2)
    X.vH=X.vH[,-1]
    
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

df  <- read.dta("sipp1991.dta")
#df = read.csv("Voting.csv")

spec=1 #spec in (1-3)
#quintile=0 #quintile in (1-5). 0 means all quintiles
quintile <- 0
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
# p number of covariables (9) including the intercept (1) and treatment (1) == 11
args = list(T = T)


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
c=0.5
alpha=0.1
tol=1e-6
set.seed(1) # for sample splitting
type = "ATE"
alpha_estimator=1
gamma_estimator=1
bias=0 #0 for unbiased or 1 for biased
L = 5 # Number of Folds
#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn, 4 2-layer nn (2-layer nn only works for ATE)
# results<-rrr(args,type,Y,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias)
get_D <- function(args,type,Y,X,m,rho_hat,b){
  #n=nrow(X)
  #p=length(b(T[1],X[1,]))
  n = get.n(args,type,X)
  p = get.p(args,type,X,b)
  
  df=matrix(0,p,n)
  if (type== "ATE"){
    T = args$T
    for (i in 1:n){
      df[,i]=b(T[i],X[i,])*as.vector(rho_hat %*% b(T[i],X[i,])) -m.ATE(Y[i],T[i],X[i,],b)
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
  folds <- split(sample(n, n,replace=FALSE), as.factor(1:L)) # sampling with replacement all the folds with the same n of observations
  
  Psi_tilde=numeric(0) #Initialize Psi_tilde as 0. This will become a running list of
  l = 1
  #Iterating through I_1,,,,I_n
  for (l in 1:L){
    #Breaks into two sets. variables with .l are not in X_i used for training. variables with .nl are used for calculating theta
    Y.l=Y[folds[[l]]] # data to generate the parameter
    Y.nl=Y[-folds[[l]]] # data use to train the model alpha and gamma
    
    X.l=X[folds[[l]],]
    X.nl=X[-folds[[l]],]
    
    
    
    if (type=="ATE"){
      
      T.l=T[folds[[l]]]
      T.nl=T[-folds[[l]]]
      
      n.l=length(T.l)
      n.nl=length(T.nl)
      
      stage1.args = list(T = T.nl)
      
      
    }
    
    # get stage 1 (on nl)
    # stage1_estimators<-get_stage1(stage1.args, type, Y.nl,X.nl,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator)
    
    get_stage1<-function(stage1.args,type,Y.nl,X.nl,p0,D_LB,D_add,max_iter,b,alpha_estimator = 1,gamma_estimator = 1){
      ###Inputs
      ###stage1.args: list of differing arguments between ATE and PD
      ###   ATE -> stage1.args = list(T)
      ###   PD -> stage1.args = list(X.nl.up,X.nl.down,delta)
      ###Y.nl = outcome variable
      ###X.nl = input variables
      ###type = ATE or PD
      ###p0: initial parameter for p
      ###D_LB: Lower bound on D
      ###D_add: value added to D
      ###max_iter: maX.nlimum iterations of chosen ML algorithm
      ###b: dictionarY.nl 
      ###alpha_estimator: numerical value indicating which model to use for the alpha estimator
      ###0 dantzig, 1 lasso
      ###gamma_estimator: numerical value indicating which model to use for the gamma estimator
      ###0 dantzig, 1 lasso, 2 rf, 3 nn
      ###Output
      ###alpha_hat and gamma_hat estimators
      
      #For ATE
      #p=length(b(T[1],X.nl[1,]))
      #n=length(T)
      
      #For PD
      #p = dim(X.nl)[2]
      #n = dim(X.nl)[1]
      
      p = get.p(stage1.args,type,X.nl,b)
      n = get.n(stage1.args,type,X.nl)
      # MNG<-get_MNG(stage1.args,type,Y.nl,X.nl,b)
      
      # Getting the M, N, G matrices
      T1 = stage1.args$T
      
      p=length(b(T1[1],X.nl[1,]))
      n.nl=length(T1)
      
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
      # B=MNG[[4]]
      
      ###########
      # alpha hat
      ###########
      if(alpha_estimator==0){ # dantzig
        
        rho_hat=RMD_stable(stage1.args,type,Y.nl,X.nl,p0,D_LB,D_add,max_iter,b,1,0)
        alpha_hat<-function(d,z){
          return(b(d,z)%*%rho_hat)
        }
        
      } else if(alpha_estimator==1){ # lasso
        
        rho_hat=RMD_stable(stage1.args,type,Y.nl,X.nl,p0,D_LB,D_add,max_iter,b,1,1)
        
        
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
          
          p = get.p(stage1.args,type,X.nl,b)
          n = get.n(stage1.args,type,X.nl)
          
          # low-dimensional moments
          X0=X[,1:p0]
          if (type == "ATE"){
            args0 = list(T = args$T)
          } else {
            args0 = list(X.up = args$X.up[,1:p0], X.down = args$X.down[,1:p0], delta = args$delta)
          }
          MNG0<-get_MNG(args0, type, Y.nl,X0,b)
          M_hat0=MNG0[[1]]
          N_hat0=MNG0[[2]]
          G_hat0=MNG0[[3]]
          
          # initial estimate
          rho_hat0=solve(G_hat0,M_hat0)
          rho_hat=c(rho_hat0,rep(0,p-ncol(G_hat0)))
          beta_hat0=solve(G_hat0,N_hat0)
          beta_hat=c(beta_hat0,rep(0,p-ncol(G_hat0)))
          
          # moments
          MNG<-get_MNG(stage1.args,type,Y.nl,X.nl,b)
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
              D_hat_rho=get_D(stage1.args,type,Y.nl,X.nl,m,rho_hat_old,b)
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
        
        
        
        
        alpha_hat<-function(d,z){
          return(b(d,z)%*%rho_hat)
        }
        
      }
      
      ###########
      # gamma hat
      ###########
      if(gamma_estimator==0){ # dantzig
        
        beta_hat=RMD_stable(stage1.args,type,Y.nl,X.nl,p0,D_LB,D_add,max_iter,b,0,0)
        gamma_hat<-function(d,z){
          return(b(d,z)%*%beta_hat)
        }
        
      } else if(gamma_estimator==1){} else if(gamma_estimator==2){ # random forest
        
        forest<- do.call(randomForest, append(list(X.nl=B,Y.nl=Y.nl), arg_Forest))
        gamma_hat<-function(d,z){
          return(predict(forest,newdata=b(d,z), type="response"))
        }
        
      } else if(gamma_estimator==3){ # neural net
        
        # scale down, de-mean, run NN, scale up, remean so that NN works well
        maX.nls_B <- apply(B, 2, maX.nl)
        mins_B <- apply(B, 2, min)
        
        maX.nls_Y.nl<-maX.nl(Y.nl)
        mins_Y.nl<-min(Y.nl)
        
        # hack to ensure that constant covariates do not become NA in the scaling
        const=maX.nls_B==mins_B
        keep=(1-const)*1:length(const)
        
        NN_B<-B
        NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maX.nls_B[keep] - mins_B[keep])
        
        NN_Y.nl<-scale(Y.nl, center = mins_Y.nl, scale = maX.nls_Y.nl - mins_Y.nl)
        
        nn<- do.call(nnet, append(list(X.nl=NN_B,Y.nl=NN_Y.nl), arg_Nnet))
        gamma_hat<-function(d,z){
          
          test<-t(as.vector(b2(d,z)))
          NN_b<-test
          NN_b[,keep]<-scale(t(NN_b[,keep]), 
                             center = mins_B[keep], 
                             scale = maX.nls_B[keep] - mins_B[keep])
          
          NN_Y.nl_hat<-predict(nn,newdata=NN_b)
          Y.nl_hat=NN_Y.nl_hat*(maX.nls_Y.nl-mins_Y.nl)+mins_Y.nl
          
          return(Y.nl_hat)
        }
        
      } else if(gamma_estimator==4){ # 2 laY.nler NN (keras)
        
        
        # scale down, de-mean, run NN, scale up, remean so that NN works well
        # hack to ensure that constant covariates do not become NA in the scaling
        maX.nls_B <- apply(B, 2, maX.nl)
        mins_B <- apply(B, 2, min)
        maX.nls_Y.nl<-maX.nl(Y.nl)
        mins_Y.nl<-min(Y.nl)
        const=maX.nls_B==mins_B
        keep=(1-const)*1:length(const)
        NN_B<-B
        NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maX.nls_B[keep] - mins_B[keep])
        NN_Y.nl<-scale(Y.nl, center = mins_Y.nl, scale = maX.nls_Y.nl - mins_Y.nl)
        
        # choose architecture
        build_model <- function() {
          model <- keras_model_sequential() %>% 
            laY.nler_dense(units = 8, activation = "relu", 
                        input_shape = dim(NN_B)[[2]]) %>% 
            laY.nler_dense(units = 8, activation = "relu") %>% 
            laY.nler_dense(units = 1) 
          
          model %>% compile(
            optimizer = "rmsprop", 
            loss = "mse", 
            metrics = c("mae")
          )
        }
        
        # use package
        model <- build_model()
        num_epochs <- 100
        model %>% fit(NN_B, NN_Y.nl,
                      epochs = num_epochs, batch_size = 1, verbose = 0)
        
        gamma_hat<-function(d,z){
          
          NN_b<-t(as.vector(b(d,z))) # test point
          NN_b[,keep]<-scale(t(NN_b[,keep]), 
                             center = mins_B[keep], 
                             scale = maX.nls_B[keep] - mins_B[keep])
          
          # 2 laY.nler NN (keras)
          NN_Y.nl_hat<-model %>% predict(NN_b, verbose = 0)
          
          Y.nl_hat=NN_Y.nl_hat*(maX.nls_Y.nl-mins_Y.nl)+mins_Y.nl
          
          return(Y.nl_hat)
        }
        
      }
      
      
      
      return(list(alpha_hat,gamma_hat))
      
    }
    
    
    alpha_hat=stage1_estimators[[1]]
    gamma_hat=stage1_estimators[[2]]
    
    print(paste0('fold: ',l))
    
    #get stage 2 (on l)
    #psi_star
    Psi_tilde.l=rep(0,n.l)
    for (i in 1:n.l){
      if (type=="ATE"){psi.args = list(T = T.l[i])} 
      else {psi.args = list(X.up = X.up.l[i,],X.down = X.down.l[i,],delta = delta)}
      if(bias){ #plug-in
        Psi_tilde.l[i]=psi_tilde_bias(psi.args, type, Y.l[i],X.l[i,],alpha_hat,gamma_hat) # without subtracting theta_hat
      }else{ #DML
        Psi_tilde.l[i]=psi_tilde(psi.args, type, Y.l[i],X.l[i,],alpha_hat,gamma_hat) # without subtracting theta_hat
      }
    }
    
    Psi_tilde=c(Psi_tilde,Psi_tilde.l)
    
    
  }
  
  #point estimation
  ate=mean(Psi_tilde)
  
  #influences
  Psi=Psi_tilde-ate
  
  var=mean(Psi^2)
  se=sqrt(var/n)
  #Returns ATE and SE
  if (type == "ATE"){
    out<-c(table(T)[[2]],table(T)[[1]],ate,se)
  } else {
    out<-c(n,ate,se)
  }
  
  return(out)
}









