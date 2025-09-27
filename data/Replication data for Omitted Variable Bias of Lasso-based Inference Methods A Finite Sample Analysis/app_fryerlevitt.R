###################################################################################################
# Empirical study: Fryer and Levitt (2013)
# Authors: Kaspar Wuthrich and Ying Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
###################################################################################################

rm(list = ls())

setwd("/Users/kasparwuthrich/Dropbox/research/Kaspar_Ying_research/submission/ReStat/final submission/replication_package_final")

library(haven)
library(xtable)
library(hdm)
library(glmnet)
library(sandwich)

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
data_test_scores <- read_dta("data_test_scores.dta")
attach(data_test_scores)

### Define specification
y <- stand_fullscale_iq_7years
d <- black 
X <- cbind(dage7y_1, 
          dage7y_2,
          dage7y_3, 
          dage7y_4, 
          dage7y_5, 
          dage7y_6, 
          dage7y_7, 
          dage7y_8, 
          dage7y_9, 
          dage7y_10, 
          dage7y_11, 
          dage7y_12, 
          dage7y_13, 
          dage7y_14, 
          dage7y_15, 
          dage7y_16, 
          dage7y_17, 
          dage7y_18, 
          female, 
          dad_hs_dropout, 
          dad_hs_grad,
          dad_some_college,
          dad_college_plus,
          dad_no_occupation,
          dad_professional,
          dad_non_professional,
          mom_hs_dropout,
          mom_hs_grad,
          mom_some_college,
          mom_college_plus,
          mom_no_occupation,
          mom_professional,
          mom_non_professional,
          inc_less_500,
          inc_500_1000,
          inc_1000_1500,
          inc_1500_2000,
          inc_2000_2500,
          inc_2500_plus,
          siblings_0,
          siblings_1,
          siblings_2,
          siblings_3,
          siblings_4,
          siblings_5,
          siblings_6_plus,
          both_bio_parents,
          age_mom,
          age_mom_2,
          age_mom_3,
          age_mom_4,
          age_mom_5,
          miss_age_mom,
          mother_indifferent,
          mother_accepting,
          mother_attentive,
          mother_over_caring,
          mother_other,
          miss_parental_score,
          w_less_1500,
          w_1500_2500,
          w_2500_3500,
          w_3500_more,
          weeks_premature_0,
          weeks_premature_1,
          weeks_premature_2,
          weeks_premature_3,
          weeks_premature_4,
          weeks_premature_5,
          weeks_premature_6,
          weeks_premature_7,
          weeks_premature_8,
          weeks_premature_9,
          weeks_premature_10,
          weeks_premature_11,
          singleton,
          twin,
          high_order_multiple)

# Test regression to compare to Stata

summary(lm(y ~ d + X))
summary(lm(y ~ d))

### Estimates

# OLS with and without covariates

ols     <- reg_robust(y,d,X)
ols_wc  <- reg_robust(y,d,X,constonly=1)

# PDL BCCH

pdl_bcch <- rlassoEffect(X,y,d,method="double selection")
res_bcch <- summary(pdl_bcch)$coef[,1:2]

# PDL BCCH x0.5

pdl_bcch_05 <- rlassoEffect(X,y,d,method="double selection",penalty=list(c = 1.1*0.5))
res_bcch_05 <- summary(pdl_bcch_05)$coef[,1:2]

# PDL BCCH x1.5

pdl_bcch_15 <- rlassoEffect(X,y,d,method="double selection",penalty=list(c = 1.1*1.5))
res_bcch_15 <- summary(pdl_bcch_15)$coef[,1:2]

###################################################################################################
# Table
###################################################################################################

alphas <- as.matrix(c(res_bcch[1],res_bcch_05[1],res_bcch_15[1],ols$alpha_hat,ols_wc$alpha_hat))
ses    <- as.matrix(c(res_bcch[2],res_bcch_05[2],res_bcch_15[2],ols$se_hat,ols_wc$se_hat))

rownames(alphas) <- c("PDL","PDL05","PDL15","OLS with covariates","OLS without covariates")
xtable(cbind(alphas,ses),digits=4)

###################################################################################################
# Simulations
###################################################################################################

nreps <- 1000

ns_vec <- c(200,400,800,1600)

res_mean <- res_sd <- matrix(NA,length(ns_vec),4)

a <- Sys.time()

for (i in 1:length(ns_vec)){
  
  ns  <- ns_vec[i]
  print(ns)
  
  obj <- sim_wrapper(y,d,X,ns,nreps)

  res_mean[i,]  <- c(ns,obj$res_mean)
  res_sd[i,]    <- c(ns,obj$res_sd)

}  

Sys.time() - a


filename <- paste("app_fryerlevitt",format(Sys.time(),"%Y-%m-%d_%T"),".RData",sep="")
save.image(file=filename)

###################################################################################################
# Figures
###################################################################################################

### Preparation

alphas_bcch   <- c(res_bcch[1],res_bcch_05[1],res_bcch_15[1])

fig_bias      <- res_mean[,-1]-rbind(alphas_bcch,alphas_bcch,alphas_bcch,alphas_bcch)
fig_bias_sd   <- (res_mean[,-1]-rbind(alphas_bcch,alphas_bcch,alphas_bcch,alphas_bcch))/res_sd[,-1]

### Figures

pdf("fryerlevitt_bias.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(ns_vec),c(-0.15,0.05), ylab="Bias", xlab="n", main="Bias", type="n")
lines(ns_vec,fig_bias[,1],type="b",pch=1,lwd=3)
lines(ns_vec,fig_bias[,2],type="b",pch=2,lwd=3)
lines(ns_vec,fig_bias[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression("0.5"*lambda[BCCH]), expression("1.5"*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()


pdf("fryerlevitt_bias_sd.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(ns_vec),c(-1.5,0.5), ylab="Bias/Std", xlab="n", main="Bias/standard deviation", type="n")
lines(ns_vec,fig_bias_sd[,1],type="b",pch=1,lwd=3)
lines(ns_vec,fig_bias_sd[,2],type="b",pch=2,lwd=3)
lines(ns_vec,fig_bias_sd[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression("0.5"*lambda[BCCH]), expression("1.5"*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

