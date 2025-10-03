######################################################################################################
# Code for reproducing the simulation results in the paper and appendix
# Authors: Kaspar Wuthrich and Ying Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
######################################################################################################

rm(list = ls())

setwd("C:/Users/187697/Documents/GitHub/master_thesis")

library(hdm)
library(glmnet)
library(sandwich)
library(foreach)
library(doParallel)
library(DoubleML)
library(mlr3)
library(mlr3learners)


# ncores <- detectCores()-1
ncores <- 10
cl     <- makeCluster(ncores)
registerDoParallel(cl)


######################################################################################################
# Functions
######################################################################################################

### Generic functions

source("code/Rfunctions.R")
source("code/adml_ovb.R")

### Simulation wrappers

# Simulations main DGPs (assumes alpha=0 and n/2>p)
sim_wrapper_main <- function(n,p,k,R2s,nreps){
  
  q_95 <- qnorm(0.95)
  res_bias <- res_cov <- res_sd <- res_leng <- res_nsel <- res_nsel_rel <- matrix(NA,length(R2s),8)
  nna <- rep(NA,length(R2s))
  
  for (j in 1:length(R2s)){
    
    temp_cov  <- temp_alpha <- temp_leng <- temp_nsel  <- temp_nsel_rel <- matrix(NA,nreps,8)
    
    R2 <- R2s[j]
    
    for (r in 1:nreps){
      
      sim_data <- gen_design(n,p,k,R2,R2,alpha0=0,dgp="homosc")
      
      y <- sim_data$y
      d <- sim_data$d
      X <- sim_data$X
      
      pdl_bcch            <- rlassoEffect(X,y,d,method="double selection")
      res_bcch            <- summary(pdl_bcch)$coef[,1:2]
      temp_alpha[r,1]     <- res_bcch[1]
      temp_cov[r,1]       <- (abs(res_bcch[1]/res_bcch[2])<=q_95)
      temp_leng[r,1]      <- res_bcch[2] 
      temp_nsel[r,1]      <- sum(pdl_bcch$selection.index)
      temp_nsel_rel[r,1]  <- sum(pdl_bcch$selection.index[1:k])
      
      pdl_min             <- pdl_cv(y,d,X,cv="min")
      temp_alpha[r,2]     <- pdl_min$alpha_hat
      temp_cov[r,2]       <- (abs(pdl_min$alpha_hat/pdl_min$se_hat)<=q_95)
      temp_leng[r,2]      <- pdl_min$se_hat
      temp_nsel[r,2]      <- sum(pdl_min$sel_index)
      temp_nsel_rel[r,2]  <- sum(pdl_min$sel_index[1:k])
      
      pdl_1se             <- pdl_cv(y,d,X,cv="1se")
      temp_alpha[r,3]     <- pdl_1se$alpha_hat
      temp_cov[r,3]       <- (abs(pdl_1se$alpha_hat/pdl_1se$se_hat)<=q_95)
      temp_leng[r,3]      <- pdl_1se$se_hat
      temp_nsel[r,3]      <- sum(pdl_1se$sel_index)
      temp_nsel_rel[r,3]  <- sum(pdl_1se$sel_index[1:k])
      
      
      
      lassodml_1se             <- dml_lasso_cf_once(y,d,X,cv_rule="1se")
      temp_alpha[r,3]     <- lassodml_1se$alpha_hat
      temp_cov[r,3]       <- (abs(lassodml_1se$alpha_hat/lassodml_1se$se_hat)<=q_95)
      temp_leng[r,3]      <- lassodml_1se$se_hat
      temp_nsel[r,3]      <- 0
      temp_nsel_rel[r,3]  <- 0
      
      
      obj_hck           <- try(reg_hck(y,d,X))[1:2]
      if (is.numeric(obj_hck)==TRUE){
        temp_alpha[r,4]   <- obj_hck[1]
        temp_cov[r,4]     <- (abs(obj_hck[1]/obj_hck[2])<=q_95)
        temp_leng[r,4]    <- obj_hck[2]
      } else {
        temp_alpha[r,4]   <- NA
        temp_cov[r,4]     <- NA
        temp_leng[r,4]    <- NA
      }
      
      pdl_bcch05          <- rlassoEffect(X,y,d,method="double selection",penalty=list(c = (1.1*0.5)))
      res_bcch05          <- summary(pdl_bcch05)$coef[,1:2]
      temp_alpha[r,5]     <- res_bcch05[1]
      temp_cov[r,5]       <- (abs(res_bcch05[1]/res_bcch05[2])<=q_95)
      temp_leng[r,5]      <- res_bcch05[2]
      temp_nsel[r,5]      <- sum(pdl_bcch05$selection.index)
      temp_nsel_rel[r,5]  <- sum(pdl_bcch05$selection.index[1:k])
      
      pdl_bcch15          <- rlassoEffect(X,y,d,method="double selection",penalty=list(c = (1.1*1.5)))
      res_bcch15          <- summary(pdl_bcch15)$coef[,1:2]
      temp_alpha[r,6]     <- res_bcch15[1]
      temp_cov[r,6]       <- (abs(res_bcch15[1]/res_bcch15[2])<=q_95)
      temp_leng[r,6]      <- res_bcch15[2]
      temp_nsel[r,6]      <- sum(pdl_bcch15$selection.index)
      temp_nsel_rel[r,6]  <- sum(pdl_bcch15$selection.index[1:k])
      
      db_min           <- db_cv(y,d,X,cv="min")
      temp_alpha[r,7]  <- db_min$alpha_hat
      temp_cov[r,7]    <- (abs(db_min$alpha_hat/db_min$se_hat)<=q_95)
      temp_leng[r,7]   <- db_min$se_hat      
      
      db_1se           <- db_cv(y,d,X,cv="1se")
      temp_alpha[r,8]  <- db_1se$alpha_hat
      temp_cov[r,8]    <- (abs(db_1se$alpha_hat/db_1se$se_hat)<=q_95)
      temp_leng[r,8]   <- db_1se$se_hat
      
      
    }
    
    set.seed(r)
    res_bias[j,]      <- colMeans(temp_alpha,na.rm=TRUE)
    res_cov[j,]       <- colMeans(temp_cov,na.rm=TRUE)
    res_sd[j,]        <- apply(temp_alpha,MARGIN=2,FUN=sd,na.rm=TRUE)
    res_leng[j,]      <- 2*q_95*colMeans(temp_leng,na.rm=TRUE)
    nna[j]            <- mean(is.na(temp_alpha[,4]))
    res_nsel[j,]      <- colMeans(temp_nsel,na.rm=TRUE)
    res_nsel_rel[j,]  <- colMeans(temp_nsel_rel,na.rm=TRUE)
    
  }
  
  return(list(res_bias=res_bias,res_cov=res_cov,res_sd=res_sd,res_leng=res_leng,res_nsel=res_nsel,res_nsel_rel=res_nsel_rel,nna=nna))
  
}
# Simulations additional DGPs
sim_wrapper_app <- function(n,p,k,R2s,alpha0,dgp,nreps){
  
  q_95 <- qnorm(0.95)
  res_bias <- res_cov <- res_sd <- res_leng <- matrix(NA,length(R2s),4)
  nna <- rep(NA,length(R2s))
  
  for (j in 1:length(R2s)){
    
    temp_cov <- temp_bias <- temp_alpha <- temp_leng <- matrix(NA,nreps,4)
    
    R2 <- R2s[j]
    for (r in 1:nreps){
      sim_data <- gen_design(n,p,k,R2,R2,alpha0=alpha0,dgp=dgp)
      
      y <- sim_data$y
      d <- sim_data$d
      X <- sim_data$X
      
      pdl_bcch          <- rlassoEffect(X,y,d,method="double selection")
      res_bcch          <- summary(pdl_bcch)$coef[,1:2]
      temp_alpha[r,1]   <- res_bcch[1]
      temp_cov[r,1]     <- (abs((res_bcch[1]-alpha0)/res_bcch[2])<=q_95)
      temp_bias[r,1]    <- res_bcch[1]-alpha0
      temp_leng[r,1]    <- res_bcch[2] 
      
      pdl_min           <- pdl_cv(y,d,X,cv="min")
      temp_alpha[r,2]   <- pdl_min$alpha_hat
      temp_cov[r,2]     <- (abs((pdl_min$alpha_hat-alpha0)/pdl_min$se_hat)<=q_95)
      temp_bias[r,2]    <- pdl_min$alpha_hat-alpha0
      temp_leng[r,2]    <- pdl_min$se_hat
      
      pdl_1se           <- pdl_cv(y,d,X,cv="1se")
      temp_alpha[r,3]   <- pdl_1se$alpha_hat
      temp_cov[r,3]     <- (abs((pdl_1se$alpha_hat-alpha0)/pdl_1se$se_hat)<=q_95)
      temp_bias[r,3]    <- pdl_1se$alpha_hat-alpha0
      temp_leng[r,3]    <- pdl_1se$se_hat
      
      if (n/2>p){
        obj_hck           <- try(reg_hck(y,d,X))[1:2]
        if (is.numeric(obj_hck)==TRUE){
          temp_alpha[r,4]   <- obj_hck[1]
          temp_cov[r,4]     <- (abs((obj_hck[1]-alpha0)/obj_hck[2])<=q_95)
          temp_bias[r,4]    <- obj_hck[1]-alpha0
          temp_leng[r,4]    <- obj_hck[2]
        } else {
          temp_alpha[r,4]   <- NA
          temp_cov[r,4]     <- NA
          temp_bias[r,4]    <- NA
          temp_leng[r,4]    <- NA
        }
      } else {
        temp_alpha[r,4]   <- NA
        temp_cov[r,4]     <- NA
        temp_bias[r,4]    <- NA
        temp_leng[r,4]    <- NA
      }
    }
    
    res_bias[j,]  <- colMeans(temp_bias,na.rm=TRUE)
    res_cov[j,]   <- colMeans(temp_cov,na.rm=TRUE)
    res_sd[j,]    <- apply(temp_alpha,MARGIN=2,FUN=sd,na.rm=TRUE)
    res_leng[j,]  <- 2*q_95*colMeans(temp_leng,na.rm=TRUE)
    nna[j]        <- mean(is.na(temp_bias[,4]))
    
  }
  
  return(list(res_bias=res_bias,res_cov=res_cov,res_sd=res_sd,res_leng=res_leng,nna=nna))
  
}  
# Simulations main DGPs (assumes alpha=0 and n/2>p) new one
restructure_results <- function(res_list) {
  n_outer <- length(res_list)          # 20
  n_mid   <- length(res_list[[1]])     # 7
  n_inner <- length(res_list[[1]][[1]])# 5
  
  # names
  mid_names   <- names(res_list[[1]])
  inner_names <- names(res_list[[1]][[1]])
  
  out <- vector("list", n_mid)
  names(out) <- mid_names
  
  for (m in seq_len(n_mid)) {
    out[[m]] <- vector("list", n_inner)
    names(out[[m]]) <- inner_names
    
    for (i in seq_len(n_inner)) {
      # extract the i-th element of the m-th slot across all outer reps
      pieces <- lapply(res_list, function(x) {
        # ensure row binding works — coerce to row if vector
        as.data.frame(x[[m]][[i]])
      })
      out[[m]][[i]] <- do.call(rbind, pieces)
    }
  }
  
  out
}
summarize_rtlist <- function(rtlist, q_95 = qnorm(0.95), alpha= 0) {
  n_mid <- length(rtlist)   # should be 7
  ncols = ncol(rtlist[[1]][["temp_alpha"]])
  res_bias     <- matrix(NA, n_mid, ncols)
  res_cov      <- matrix(NA, n_mid, ncols)
  res_sd       <- matrix(NA, n_mid, ncols)
  res_leng     <- matrix(NA, n_mid, ncols)
  res_nsel     <- matrix(NA, n_mid, ncols)
  res_nsel_rel <- matrix(NA, n_mid, ncols)
  nna          <- rep(NA, n_mid)
  
  for (j in seq_len(n_mid)) {
    temp_alpha    <- as.matrix(rtlist[[j]][["temp_alpha"]])-alpha
    temp_leng     <- as.matrix(rtlist[[j]][["temp_leng"]])
    tval <-  qnorm(1 - 0.05/2)
    # confidence intervals cell by cell
    ci_low  <- as.matrix(rtlist[[j]][["temp_alpha"]]) - tval * temp_leng
    ci_high <- as.matrix(rtlist[[j]][["temp_alpha"]]) + tval * temp_leng
    temp_cov <- (ci_low <= alpha & ci_high >= alpha)
    
    temp_nsel     <- as.matrix(rtlist[[j]][["temp_nsel"]])
    temp_nsel_rel <- as.matrix(rtlist[[j]][["temp_nsel_rel"]])
    
    res_bias[j, ]     <- colMeans(abs(temp_alpha), na.rm = TRUE)
    res_cov[j, ]      <- colMeans(temp_cov, na.rm = TRUE)
    res_sd[j, ]       <- apply(temp_alpha,MARGIN=2,FUN=sd,na.rm=TRUE)
    res_leng[j, ]     <- 2*q_95*colMeans(temp_leng,na.rm=TRUE)
    nna[j]            <- mean(is.na(temp_alpha[, 4]))
    res_nsel[j, ]     <- colMeans(temp_nsel, na.rm = TRUE)
    res_nsel_rel[j, ] <- colMeans(temp_nsel_rel, na.rm = TRUE)
  }
  
  list(
    res_bias = res_bias,
    res_cov = res_cov,
    res_sd = res_sd,
    res_leng = res_leng,
    nna = nna,
    res_nsel = res_nsel,
    res_nsel_rel = res_nsel_rel
  )
}
sim_wrapper_main_aux <- function(n,p,k,R2s, seedval = 1, nreps=1){
  
  
  alphaval <- 0.05
  q_95 <- qnorm(1 - alphaval/2)
  r2info <- list()
  for (j in 1:length(R2s)){
    
    temp_cov  <- temp_alpha <- temp_leng <- temp_nsel  <- temp_nsel_rel <- matrix(NA,1,10)
    
    R2 <- R2s[j]
    r <- 1
    # for (r in 1:nreps){
    
    sim_data <- gen_design(n,p,k,R2,R2,alpha0=4,dgp="homosc", seedval = seedval)
    
    y <- sim_data$y
    d <- sim_data$d
    X <- sim_data$X
    
    pdl_bcch            <- rlassoEffect(X,y,d,method="double selection")
    res_bcch            <- summary(pdl_bcch)$coef[,1:2]
    temp_alpha[r,1]     <- res_bcch[1]
    temp_cov[r,1]       <- (abs(res_bcch[1]/res_bcch[2])<=q_95)
    temp_leng[r,1]      <- res_bcch[2] 
    temp_nsel[r,1]      <- sum(pdl_bcch$selection.index)
    temp_nsel_rel[r,1]  <- sum(pdl_bcch$selection.index[1:k])
    
    pdl_min             <- pdl_cv(y,d,X,cv="min")
    temp_alpha[r,2]     <- pdl_min$alpha_hat
    temp_cov[r,2]       <- (abs(pdl_min$alpha_hat/pdl_min$se_hat)<=q_95)
    temp_leng[r,2]      <- pdl_min$se_hat
    temp_nsel[r,2]      <- sum(pdl_min$sel_index)
    temp_nsel_rel[r,2]  <- sum(pdl_min$sel_index[1:k])
    
    pdl_1se             <- pdl_cv(y,d,X,cv="1se")
    temp_alpha[r,3]     <- pdl_1se$alpha_hat
    temp_cov[r,3]       <- (abs(pdl_1se$alpha_hat/pdl_1se$se_hat)<=q_95)
    temp_leng[r,3]      <- pdl_1se$se_hat
    temp_nsel[r,3]      <- sum(pdl_1se$sel_index)
    temp_nsel_rel[r,3]  <- sum(pdl_1se$sel_index[1:k])
    
    obj_hck           <- try(reg_hck(y,d,X))[1:2]
    if (is.numeric(obj_hck)==TRUE){
      temp_alpha[r,4]   <- obj_hck[1]
      temp_cov[r,4]     <- (abs(obj_hck[1]/obj_hck[2])<=q_95)
      temp_leng[r,4]    <- obj_hck[2]
    } else {
      temp_alpha[r,4]   <- NA
      temp_cov[r,4]     <- NA
      temp_leng[r,4]    <- NA
    }
    
    pdl_bcch05          <- rlassoEffect(X,y,d,method="double selection",penalty=list(c = (1.1*0.5)))
    res_bcch05          <- summary(pdl_bcch05)$coef[,1:2]
    temp_alpha[r,5]     <- res_bcch05[1]
    temp_cov[r,5]       <- (abs(res_bcch05[1]/res_bcch05[2])<=q_95)
    temp_leng[r,5]      <- res_bcch05[2]
    temp_nsel[r,5]      <- sum(pdl_bcch05$selection.index)
    temp_nsel_rel[r,5]  <- sum(pdl_bcch05$selection.index[1:k])
    
    pdl_bcch15          <- rlassoEffect(X,y,d,method="double selection",penalty=list(c = (1.1*1.5)))
    res_bcch15          <- summary(pdl_bcch15)$coef[,1:2]
    temp_alpha[r,6]     <- res_bcch15[1]
    temp_cov[r,6]       <- (abs(res_bcch15[1]/res_bcch15[2])<=q_95)
    temp_leng[r,6]      <- res_bcch15[2]
    temp_nsel[r,6]      <- sum(pdl_bcch15$selection.index)
    temp_nsel_rel[r,6]  <- sum(pdl_bcch15$selection.index[1:k])
    
    db_min              <- db_cv(y,d,X,cv="min")
    temp_alpha[r,7]     <- db_min$alpha_hat
    temp_cov[r,7]       <- (abs(db_min$alpha_hat/db_min$se_hat)<=q_95)
    temp_leng[r,7]      <- db_min$se_hat      
    
    db_1se              <- db_cv(y,d,X,cv="1se")
    temp_alpha[r,8]     <- db_1se$alpha_hat
    temp_cov[r,8]       <- (abs(db_1se$alpha_hat/db_1se$se_hat)<=q_95)
    temp_leng[r,8]      <- db_1se$se_hat

    
    learner = lrn("regr.cv_glmnet", s="lambda.min")
    dml_data_sim = double_ml_data_from_matrix(X=X, y=y, d=d)
    ml_l_sim = learner$clone()
    ml_m_sim = learner$clone()
    obj_dml_plr_sim = DoubleMLPLR$new(dml_data_sim, ml_l=ml_l_sim, ml_m=ml_m_sim)
    obj_dml_plr_sim$fit()
    
    temp_alpha[r,9]     <- obj_dml_plr_sim$all_coef
    temp_cov[r,9]       <- (abs(obj_dml_plr_sim$all_coef/obj_dml_plr_sim$all_se)<=q_95)
    temp_leng[r,9]      <- obj_dml_plr_sim$all_se
    
    adml.res <- adml.lasso( sim_data )
    temp_alpha[r,10]     <- adml.res$coef
    temp_cov[r,10]       <- (abs(adml.res$coef/adml.res$se)<=q_95)
    temp_leng[r,10]      <- adml.res$se
    
    r2info[[j]] <- list(temp_alpha=temp_alpha,
                        temp_leng=temp_leng ,
                        temp_nsel=temp_nsel,
                        temp_nsel_rel=temp_nsel_rel)
  }
  return(r2info)
}


######################################################################################################



######################################################################################################
# Simulations
######################################################################################################

### Setup 

nreps <- 20
R2s       <- 1-c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)
R2s_bern  <- R2s

### Main text
a <- Sys.time()
res_list <- foreach(b = 1:nreps,
                    .packages = c("hdm", "glmnet", "sandwich", "DoubleML", "mlr3", "mlr3learners", 
                                  "foreign", "dplyr", "ggplot2", "quantreg", 
                                  "randomForest", "keras", "data.table"))  %dopar% {
  sim_wrapper_main_aux(n=500,p=200,k=5,R2s=R2s, seedval = b)
}
dgp_m1 <- summarize_rtlist( restructure_results( res_list ), alpha =4 )
Sys.time() - a


### Main text
a <- Sys.time()
res_list2 <- foreach(b = 1:nreps,
                    .packages = c("hdm", "glmnet", "sandwich", "DoubleML", "mlr3", "mlr3learners")
) %dopar% {
  sim_wrapper_main_aux(n=1000,p=200,k=5,R2s=R2s, seedval = b)
}
dgp_m2 <- summarize_rtlist( restructure_results( res_list2 ), alpha =4 )
Sys.time() - a


### Main text
a <- Sys.time()
res_list <- foreach(b = 1:nreps,
                    .packages = c("hdm", "glmnet", "sandwich", "DoubleML", "mlr3", "mlr3learners")
) %dopar% {
  sim_wrapper_main_aux(n=5000,p=200,k=5,R2s=R2s, seedval = b)
}
dgp_m3 <- summarize_rtlist( restructure_results( res_list ), alpha =4 )
Sys.time() - a


# Save objects as RDS files
saveRDS(dgp_m1, file = "dgp_m1.rds")
saveRDS(dgp_m2, file = "dgp_m2.rds")
saveRDS(dgp_m3, file = "dgp_m3.rds")

######################################################################################################



################################################################################
############################## Importing Results ###############################

dgp_m1 <- readRDS("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/data/dgp_m2.rds")
dgp_m1
plot(range(R2s),c(min(dgp_m1$res_bias),max(dgp_m1$res_bias)), ylab="Bias", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_bias[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_bias[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_bias[,3],type="b",pch=3,lwd=3)
lines(R2s,dgp_m1$res_bias[,9],type="b",pch=4,lwd=3)
lines(R2s,dgp_m1$res_bias[,10],type="b",pch=5,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)


dgp_m1 <- readRDS("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/data/dgp_m1.rds")
plot(range(R2s),c(0.91,0.96), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_cov[,3],type="b",pch=3,lwd=3)
# lines(R2s,dgp_m1$res_cov[,9],type="b",pch=4,lwd=3)
lines(R2s,dgp_m1$res_cov[,10],type="b",pch=5,lwd=3)


dgp_m1 <- readRDS("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/data/dgp_m2.rds")
plot(range(R2s),c(0.91,0.97), ylab="Coverage", xlab=expression(R^{2}), main="n=5000, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_cov[,3],type="b",pch=3,lwd=3)
# lines(R2s,dgp_m1$res_cov[,9],type="b",pch=4,lwd=3)
lines(R2s,dgp_m1$res_cov[,10],type="b",pch=5,lwd=3)



dgp_m1 <- readRDS("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/data/dgp_m1.rds")
dgp_m2 <- readRDS("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/data/dgp_m2.rds")
dgp_m3 <- readRDS("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/data/dgp_m3.rds")
plot(range(R2s),c(0.91,0.97), ylab="Coverage", xlab=expression(R^{2}), main="n=5000, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_cov[,10],type="b",pch=5,lwd=3)
lines(R2s,dgp_m2$res_cov[,10],type="b",pch=5,lwd=3)
lines(R2s,dgp_m3$res_cov[,10],type="b",pch=5,lwd=3)



png("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/output/coverage_plot.png", width=800, height=600, res=120)
plot(range(R2s), c(0.91,1),
     ylab="Coverage",
     xlab=expression(R^{2}),
     main="p=200, k=5",
     type="n")

# First 3 lines (blue)
lines(R2s, dgp_m1$res_cov[,1],  type="b", pch=1, lwd=3, col="blue")
lines(R2s, dgp_m2$res_cov[,1],  type="b", pch=2, lwd=3, col="blue")
lines(R2s, dgp_m3$res_cov[,1],  type="b", pch=3, lwd=3, col="blue")

# Last 3 lines (red)
lines(R2s, dgp_m1$res_cov[,10], type="b", pch=1, lwd=3, col="red")
lines(R2s, dgp_m2$res_cov[,10], type="b", pch=2, lwd=3, col="red")
lines(R2s, dgp_m3$res_cov[,10], type="b", pch=3, lwd=3, col="red")

lines(R2s, dgp_m1$res_cov[,2], type="b", pch=1, lwd=3, col="green")
lines(R2s, dgp_m2$res_cov[,2], type="b", pch=2, lwd=3, col="green")
lines(R2s, dgp_m3$res_cov[,2], type="b", pch=3, lwd=3, col="green")

# Legend for pch
legend("topright",
       legend=c("n=500", "n=1000", "n=5000"),
       pch=c(1, 2, 3),
       bty="n")
legend("topleft",
       legend=c("OLS", "DML Lasso", "Double Lasso"),
       col=c("blue", "red", "green"),
       lwd=3,
       bty="n")
dev.off()








png("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/output/bias_plot.png", width=800, height=600, res=120)
plot(range(R2s), c(0.0,0.05),
     ylab="Bias",
     xlab=expression(R^{2}),
     main="p=200, k=5",
     type="n")

# First 3 lines (blue)
lines(R2s, dgp_m1$res_bias[,1],  type="b", pch=1, lwd=3, col="blue")
lines(R2s, dgp_m2$res_bias[,1],  type="b", pch=2, lwd=3, col="blue")
lines(R2s, dgp_m3$res_bias[,1],  type="b", pch=3, lwd=3, col="blue")

# Last 3 lines (red)
lines(R2s, dgp_m1$res_bias[,10], type="b", pch=1, lwd=3, col="red")
lines(R2s, dgp_m2$res_bias[,10], type="b", pch=2, lwd=3, col="red")
lines(R2s, dgp_m3$res_bias[,10], type="b", pch=3, lwd=3, col="red")

# Last 3 lines (red)
lines(R2s, dgp_m1$res_bias[,9], type="b", pch=1, lwd=3, col="black")
lines(R2s, dgp_m2$res_bias[,9], type="b", pch=2, lwd=3, col="black")
lines(R2s, dgp_m3$res_bias[,9], type="b", pch=3, lwd=3, col="black")

lines(R2s, dgp_m1$res_bias[,2], type="b", pch=1, lwd=3, col="green")
lines(R2s, dgp_m2$res_bias[,2], type="b", pch=2, lwd=3, col="green")
lines(R2s, dgp_m3$res_bias[,2], type="b", pch=3, lwd=3, col="green")

# Legend for pch
legend("topright",
       legend=c("n=500", "n=1000", "n=5000"),
       pch=c(1, 2, 3),
       bty="n")
legend("topleft",
       legend=c("OLS",  "DML", "ADML",  "Double Lasso"),
       col=c("blue", "black", "red", "green"),
       lwd=3,
       bty="n")
dev.off()





png("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/output/biasse_plot.png", width=800, height=600, res=120)
plot(range(R2s), c(0.77,0.86),
     ylab="Bias/SE",
     xlab=expression(R^{2}),
     main="p=200, k=5",
     type="n")

# First 3 lines (blue)
lines(R2s, dgp_m1$res_bias[,1] / dgp_m1$res_sd[,1],  type="b", pch=1, lwd=3, col="blue")
lines(R2s, dgp_m2$res_bias[,1] / dgp_m2$res_sd[,1],  type="b", pch=2, lwd=3, col="blue")
lines(R2s, dgp_m3$res_bias[,1] / dgp_m3$res_sd[,1],  type="b", pch=3, lwd=3, col="blue")
lines(R2s, dgp_m1$res_bias[,10] / dgp_m1$res_sd[,10], type="b", pch=1, lwd=3, col="red")
lines(R2s, dgp_m2$res_bias[,10] / dgp_m2$res_sd[,10], type="b", pch=2, lwd=3, col="red")
lines(R2s, dgp_m3$res_bias[,10] / dgp_m3$res_sd[,10], type="b", pch=3, lwd=3, col="red")
lines(R2s, dgp_m1$res_bias[,2] / dgp_m1$res_sd[,2], type="b", pch=1, lwd=3, col="green")
lines(R2s, dgp_m2$res_bias[,2] / dgp_m2$res_sd[,2], type="b", pch=2, lwd=3, col="green")
lines(R2s, dgp_m3$res_bias[,2] / dgp_m3$res_sd[,2], type="b", pch=3, lwd=3, col="green")
# Legend for pch
legend("topright",
       legend=c("n=500", "n=1000", "n=5000"),
       pch=c(1, 2, 3),
       bty="n")
legend("topleft",
       legend=c("OLS", "DML Lasso", "Double Lasso"),
       col=c("blue", "red", "green"),
       lwd=3,
       bty="n")
dev.off()

######################################################################################################


# 
# plot(range(R2s),c(min(dgp_m1$res_bias),max(dgp_m1$res_bias)), ylab="Bias", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
# lines(R2s,dgp_m1$res_bias[,1],type="b",pch=1,lwd=3)
# lines(R2s,dgp_m1$res_bias[,2],type="b",pch=2,lwd=3)
# lines(R2s,dgp_m1$res_bias[,3],type="b",pch=3,lwd=3)
# lines(R2s,dgp_m1$res_bias[,9],type="b",pch=4,lwd=3)
# lines(R2s,dgp_m1$res_bias[,10],type="b",pch=5,lwd=3)
# abline(h=0,col="darkgrey",lty=2,lwd=1)
# 
# 
# 
# res_list_shape <- restructure_results( res_list )
# 
# hist(res_list_shape[[1]][[3]]$V9)
# 
# mean(res_list_shape[[1]][[1]]$V9-4)
# 
# 
# 
# 
# plot(range(R2s),c(min(dgp_m1$res_cov),max(dgp_m1$res_cov)), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
# lines(R2s,dgp_m1$res_cov[,1],type="b",pch=1,lwd=3)
# lines(R2s,dgp_m1$res_cov[,2],type="b",pch=2,lwd=3)
# lines(R2s,dgp_m1$res_cov[,3],type="b",pch=3,lwd=3)
# lines(R2s,dgp_m1$res_cov[,9],type="b",pch=4,lwd=3)
# lines(R2s,dgp_m1$res_cov[,10],type="b",pch=5,lwd=3)



# 
# 
# dgp_m1
