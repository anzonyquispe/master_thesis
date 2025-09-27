###################################################################################################
# Empirical study: Effect of 401(k) eligibility (Belloni et al. 2017)
# Authors: Kaspar Wuthrich and Ying Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
###################################################################################################

rm(list = ls())

setwd("/Users/kasparwuthrich/Dropbox/research/Kaspar_Ying_research/submission/ReStat/final submission/replication_package_final")

library(sandwich)
library(R.matlab)
library(xtable)
library(hdm)
library(glmnet)

set.seed(12345)

###################################################################################################
# Functions
###################################################################################################

### Generic functions

source("Rfunctions.R")

### Simulation wrapper

sim_wrapper <- function(y,d,X,ns,nreps){
  
  n <- dim(X)[1]
  
  temp_res <- matrix(NA,nreps,3)
  
  for (r in 1:nreps){
    
    ind <- sample(n,ns,replace=TRUE)
    
    Xs <- X[ind,]
    ys <- y[ind]
    ds <- d[ind]
    
    pdl_bcch <- rlassoEffect(Xs,ys,ds,method="double selection")
    res_bcch <- summary(pdl_bcch)$coef[,1:2]
    
    pdl_bcch_05 <- rlassoEffect(Xs,ys,ds,method="double selection",penalty=list(c = 1.1*0.5))
    res_bcch_05 <- summary(pdl_bcch_05)$coef[,1:2]
    
    pdl_bcch_15 <- rlassoEffect(Xs,ys,ds,method="double selection",penalty=list(c = 1.1*1.5))
    res_bcch_15 <- summary(pdl_bcch_15)$coef[,1:2]
    
    temp_res[r,] <- c(res_bcch[1],res_bcch_05[1],res_bcch_15[1])
    
  }
  
  res_mean  <- colMeans(temp_res)
  res_sd    <- apply(temp_res,MARGIN=2,FUN=sd,na.rm=TRUE)
  
  return(list(res_mean=res_mean,res_sd=res_sd))
 
}

###################################################################################################
# Estimates based on the original data
###################################################################################################

### Load data

spec1 <- readMat("data_spec1.mat")
y1    <- spec1$y
d1    <- spec1$d
X1    <- spec1$Xmat

spec2 <- readMat("data_spec2.mat")
y2    <- spec2$y
d2    <- spec2$d
X2    <- spec2$Xmat

### Estimates

# OLS with and without covariates

ols1 <- reg_robust(y1,d1,X1)
ols2 <- reg_robust(y2,d2,X2)

ols1_wc  <- reg_robust(y1,d1,X1,constonly=1)
ols2_wc  <- reg_robust(y2,d2,X2,constonly=1)

# PDL BCCH

pdl_bcch1 <- rlassoEffect(X1,y1,d1,method="double selection")
res_bcch1 <- summary(pdl_bcch1)$coef[,1:2]

pdl_bcch2 <- rlassoEffect(X2,y2,d2,method="double selection")
res_bcch2 <- summary(pdl_bcch2)$coef[,1:2]

# PDL BCCH x0.5

pdl_bcch1_05 <- rlassoEffect(X1,y1,d1,method="double selection",penalty=list(c = 1.1*0.5))
res_bcch1_05 <- summary(pdl_bcch1_05)$coef[,1:2]

pdl_bcch2_05 <- rlassoEffect(X2,y2,d2,method="double selection",penalty=list(c = 1.1*0.5))
res_bcch2_05 <- summary(pdl_bcch2_05)$coef[,1:2]

# PDL BCCH x1.5

pdl_bcch1_15 <- rlassoEffect(X1,y1,d1,method="double selection",penalty=list(c = 1.1*1.5))
res_bcch1_15 <- summary(pdl_bcch1_15)$coef[,1:2]

pdl_bcch2_15 <- rlassoEffect(X2,y2,d2,method="double selection",penalty=list(c = 1.1*1.5))
res_bcch2_15 <- summary(pdl_bcch2_15)$coef[,1:2]

###################################################################################################
# Table
###################################################################################################

alphas1 <- as.matrix(c(res_bcch1[1],res_bcch1_05[1],res_bcch1_15[1],ols1$alpha_hat,ols1_wc$alpha_hat))
ses1    <- as.matrix(c(res_bcch1[2],res_bcch1_05[2],res_bcch1_15[2],ols1$se_hat,ols1_wc$se_hat))

alphas2 <- as.matrix(c(res_bcch2[1],res_bcch2_05[1],res_bcch2_15[1],ols2$alpha_hat,ols2_wc$alpha_hat))
ses2    <- as.matrix(c(res_bcch2[2],res_bcch2_05[2],res_bcch2_15[2],ols2$se_hat,ols2_wc$se_hat))

rownames(alphas1) <- c("PDL","PDL05","PDL15","OLS with covariates","OLS without covariates")
rownames(alphas2) <- c("PDL","PDL05","PDL15","OLS with covariates","OLS without covariates")

xtable(cbind(alphas1,ses1),digits=2)
xtable(cbind(alphas2,ses2),digits=2)

###################################################################################################
# Simulations
###################################################################################################

nreps <- 1000
ns_vec <- c(200,400,800,1600)

res_mean1 <- res_mean2 <- res_sd1 <- res_sd2 <- matrix(NA,length(ns_vec),4)

a<-Sys.time()

for (i in 1:length(ns_vec)){
  
  ns <- ns_vec[i]
  print(ns)
  
  obj1 <- sim_wrapper(y1,d1,X1,ns,nreps)
  obj2 <- sim_wrapper(y2,d2,X2,ns,nreps)

  res_mean1[i,]  <- c(ns,obj1$res_mean)
  res_mean2[i,]  <- c(ns,obj2$res_mean)
  
  res_sd1[i,]  <- c(ns,obj1$res_sd)
  res_sd2[i,]  <- c(ns,obj2$res_sd)
  
}  

Sys.time()-a

filename <- paste("app_401k",format(Sys.time(),"%Y-%m-%d_%T"),".RData",sep="")
save.image(file=filename)

###################################################################################################
# Figures
###################################################################################################

### Preparation

alphas_bcch1  <- c(res_bcch1[1],res_bcch1_05[1],res_bcch1_15[1])
alphas_bcch2  <- c(res_bcch2[1],res_bcch2_05[1],res_bcch2_15[1])

fig_bias1     <- res_mean1[,-1]-rbind(alphas_bcch1,alphas_bcch1,alphas_bcch1,alphas_bcch1)
fig_bias_sd1  <- (res_mean1[,-1]-rbind(alphas_bcch1,alphas_bcch1,alphas_bcch1,alphas_bcch1))/res_sd1[,-1]

fig_bias2     <- res_mean2[,-1]-rbind(alphas_bcch2,alphas_bcch2,alphas_bcch2,alphas_bcch2)
fig_bias_sd2  <- (res_mean2[,-1]-rbind(alphas_bcch2,alphas_bcch2,alphas_bcch2,alphas_bcch2))/res_sd2[,-1]

### Specification 1

pdf("401k_bias1.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(ns_vec),c(-1000,15000), ylab="Bias", xlab="n", main="Bias", type="n")
lines(ns_vec,fig_bias1[,1],type="b",pch=1,lwd=3)
lines(ns_vec,fig_bias1[,2],type="b",pch=2,lwd=3)
lines(ns_vec,fig_bias1[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression("0.5"*lambda[BCCH]), expression("1.5"*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("401k_bias_sd1.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(ns_vec),c(-0.1,1), ylab="Bias/Std", xlab="n", main="Bias/standard deviation", type="n")
lines(ns_vec,fig_bias_sd1[,1],type="b",pch=1,lwd=3)
lines(ns_vec,fig_bias_sd1[,2],type="b",pch=2,lwd=3)
lines(ns_vec,fig_bias_sd1[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression("0.5"*lambda[BCCH]), expression("1.5"*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

### Specification 2

pdf("401k_bias2.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(ns_vec),c(-4000,10000), ylab="Bias", xlab="n", main="Bias", type="n")
lines(ns_vec,fig_bias2[,1],type="b",pch=1,lwd=3)
lines(ns_vec,fig_bias2[,2],type="b",pch=2,lwd=3)
lines(ns_vec,fig_bias2[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression("0.5"*lambda[BCCH]), expression("1.5"*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("401k_bias_sd2.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(ns_vec),c(-0.4,0.8), ylab="Bias/Std", xlab="n", main="Bias/standard deviation", type="n")
lines(ns_vec,fig_bias_sd2[,1],type="b",pch=1,lwd=3)
lines(ns_vec,fig_bias_sd2[,2],type="b",pch=2,lwd=3)
lines(ns_vec,fig_bias_sd2[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression("0.5"*lambda[BCCH]), expression("1.5"*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()
