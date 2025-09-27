######################################################################################################
# Code for reproducing the simulation results in the paper and appendix
# Authors: Kaspar Wuthrich and Ying Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
######################################################################################################

rm(list = ls())

setwd("C:/Users/187697/Dropbox/data_eco/Udesa/mater_thesis/data/Replication data for Omitted Variable Bias of Lasso-based Inference Methods A Finite Sample Analysis")

library(hdm)
library(glmnet)
library(sandwich)
library(foreach)
library(doParallel)
library(DoubleML)
library(mlr3)
library(mlr3learners)
set.seed(12345)

ncores <- detectCores()-1
cl     <- makeCluster(ncores)
registerDoParallel(cl)


######################################################################################################
# Functions
######################################################################################################

### Generic functions

source("Rfunctions.R")

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
summarize_rtlist <- function(rtlist, q_95 = qnorm(0.95)) {
  n_mid <- length(rtlist)   # should be 7
  inner_names <- names(rtlist[[1]])  # e.g. temp_alpha, temp_cov, ...
  
  res_bias     <- matrix(NA, n_mid, ncol(rtlist[[1]][["temp_alpha"]]))
  res_cov      <- matrix(NA, n_mid, ncol(rtlist[[1]][["temp_cov"]]))
  res_sd       <- matrix(NA, n_mid, ncol(rtlist[[1]][["temp_alpha"]]))
  res_leng     <- matrix(NA, n_mid, ncol(rtlist[[1]][["temp_leng"]]))
  res_nsel     <- matrix(NA, n_mid, ncol(rtlist[[1]][["temp_nsel"]]))
  res_nsel_rel <- matrix(NA, n_mid, ncol(rtlist[[1]][["temp_nsel_rel"]]))
  nna          <- rep(NA, n_mid)
  
  for (j in seq_len(n_mid)) {
    temp_alpha    <- as.matrix(rtlist[[j]][["temp_alpha"]])
    temp_cov      <- as.matrix(rtlist[[j]][["temp_cov"]])
    temp_leng     <- as.matrix(rtlist[[j]][["temp_leng"]])
    temp_nsel     <- as.matrix(rtlist[[j]][["temp_nsel"]])
    temp_nsel_rel <- as.matrix(rtlist[[j]][["temp_nsel_rel"]])
    
    res_bias[j, ]     <- colMeans(temp_alpha, na.rm = TRUE)
    res_cov[j, ]      <- colMeans(temp_cov, na.rm = TRUE)
    res_sd[j, ]       <- apply(temp_alpha, 2, sd, na.rm = TRUE)
    res_leng[j, ]     <- 2 * q_95 * colMeans(temp_leng, na.rm = TRUE)
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
  
  set.seed(seedval)
  q_95 <- qnorm(0.95)
  r2info <- list()
  for (j in 1:length(R2s)){
    
    temp_cov  <- temp_alpha <- temp_leng <- temp_nsel  <- temp_nsel_rel <- matrix(NA,1,10)
    
    R2 <- R2s[j]
    r <- 1
    # for (r in 1:nreps){
      
      sim_data <- gen_design(n,p,k,R2,R2,alpha0=4,dgp="homosc")
      
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
      
      db_min           <- db_cv(y,d,X,cv="min")
      temp_alpha[r,7]  <- db_min$alpha_hat
      temp_cov[r,7]    <- (abs(db_min$alpha_hat/db_min$se_hat)<=q_95)
      temp_leng[r,7]   <- db_min$se_hat      
      
      db_1se           <- db_cv(y,d,X,cv="1se")
      temp_alpha[r,8]  <- db_1se$alpha_hat
      temp_cov[r,8]    <- (abs(db_1se$alpha_hat/db_1se$se_hat)<=q_95)
      temp_leng[r,8]   <- db_1se$se_hat
      
      
      dml_data_sim = double_ml_data_from_matrix(X=X, y=y, d=d)
      learner = lrn("regr.ranger", num.trees=100, max.depth=5, min.node.size=2)
      ml_l_bonus = learner$clone()
      ml_m_bonus = learner$clone()
      
      learner = lrn("regr.cv_glmnet", s="lambda.min")
      ml_l_sim = learner$clone()
      ml_m_sim = learner$clone()
      
      obj_dml_plr_bonus = DoubleMLPLR$new(dml_data_sim, ml_l=ml_l_bonus, ml_m=ml_m_bonus)
      obj_dml_plr_sim = DoubleMLPLR$new(dml_data_sim, ml_l=ml_l_sim, ml_m=ml_m_sim)
      obj_dml_plr_bonus$fit()
      obj_dml_plr_sim$fit()

      temp_alpha[r,9]     <- obj_dml_plr_bonus$all_coef
      temp_cov[r,9]       <- (abs(obj_dml_plr_bonus$all_coef/obj_dml_plr_bonus$all_se)<=q_95)
      temp_leng[r,9]      <- obj_dml_plr_bonus$all_se
      
      temp_alpha[r,10]    <- obj_dml_plr_sim$all_coef
      temp_cov[r  ,10]    <- (abs(obj_dml_plr_sim$all_coef/obj_dml_plr_sim$all_se)<=q_95)
      temp_leng[r ,10]    <- obj_dml_plr_sim$all_se

      
      
    # }
    
    
    
    r2info[[j]] <- list(temp_alpha=temp_alpha,
                        temp_cov=temp_cov,
                        temp_leng=temp_leng,
                        temp_nsel=temp_nsel,
                        temp_nsel_rel=temp_nsel_rel)
    
  }
  return(r2info)
  # 
  # return(list(res_bias=res_bias,res_cov=res_cov,
  #             res_sd=res_sd,res_leng=res_leng,
  #             res_nsel=res_nsel,res_nsel_rel=res_nsel_rel,nna=nna))
  # 
}


######################################################################################################



######################################################################################################
# Simulations
######################################################################################################

### Setup

nreps <- 20

R2s       <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)
R2s_bern  <- R2s

### Main text

a <- Sys.time()
res_list <- foreach(b = 1:nreps,
                    .packages = c("hdm", "glmnet", "sandwich")
) %dopar% {
  
  sim_wrapper_main_aux(n=500,p=200,k=5,R2s=R2s, seedval = b)
}
dgp_m1 <- summarize_rtlist( restructure_results( res_list ) )
# dgp_m1 <- sim_wrapper_main(n=500,p=200,k=5,R2s=R2s,nreps=nreps)
# dgp_m2 <- sim_wrapper_main(n=1000,p=200,k=5,R2s=R2s,nreps=nreps)
# dgp_m3 <- sim_wrapper_main(n=500,p=200,k=10,R2s=R2s,nreps=nreps)
Sys.time() - a

m1_bias_sd <- dgp_m1$res_bias/dgp_m1$res_sd
m2_bias_sd <- dgp_m2$res_bias/dgp_m2$res_sd
m3_bias_sd <- dgp_m3$res_bias/dgp_m3$res_sd

# Check for NAs
dgp_m1$nna
dgp_m2$nna
dgp_m3$nna

### Appendix

a <- Sys.time()
dgp_a1 <- sim_wrapper_app(n=500,p=200,k=5,R2s=R2s_bern,alpha0=0,dgp="bern",nreps=nreps)
dgp_a2 <- sim_wrapper_app(n=500,p=200,k=5,R2s=R2s,alpha0=0,dgp="tdistr",nreps=nreps)
dgp_a3 <- sim_wrapper_app(n=500,p=200,k=5,R2s=R2s,alpha0=0,dgp="heterosc",nreps=nreps)
dgp_a4 <- sim_wrapper_app(n=200,p=200,k=5,R2s=R2s,alpha0=0,dgp="homosc",nreps=nreps)
dgp_a5 <- sim_wrapper_app(n=500,p=200,k=5,R2s=R2s,alpha0=1,dgp="homosc",nreps=nreps)
dgp_a6 <- sim_wrapper_app(n=500,p=200,k=5,R2s=R2s,alpha0=-1,dgp="homosc",nreps=nreps)
Sys.time() - a

a1_bias_sd <- dgp_a1$res_bias/dgp_a1$res_sd
a2_bias_sd <- dgp_a2$res_bias/dgp_a2$res_sd
a3_bias_sd <- dgp_a3$res_bias/dgp_a3$res_sd
a4_bias_sd <- dgp_a4$res_bias/dgp_a4$res_sd
a5_bias_sd <- dgp_a5$res_bias/dgp_a5$res_sd
a6_bias_sd <- dgp_a6$res_bias/dgp_a6$res_sd

# Check for NAs
dgp_a1$nna
dgp_a2$nna
dgp_a3$nna
dgp_a4$nna
dgp_a5$nna
dgp_a6$nna

### Save RData

filename <- paste("simulations",format(Sys.time(),"%Y-%m-%d_%T"),".RData",sep="")
save.image(file=filename)

######################################################################################################
# Figures main text
######################################################################################################

### Selection

# BCCH, CV min, CV 1se

pdf("nsel_1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,40), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_nsel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_nsel[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_nsel[,3],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("right", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,40), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_nsel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_nsel[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m2$res_nsel[,3],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("right", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,60), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_nsel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_nsel[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m3$res_nsel[,3],type="b",pch=4,lwd=3)
abline(h=10,col="darkgrey",lty=2,lwd=1)
legend("right", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

# Sensitivity to regularization choice

pdf("nsel_1_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,40), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_nsel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_nsel[,5],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_nsel[,6],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_2_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,40), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_nsel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_nsel[,5],type="b",pch=2,lwd=3)
lines(R2s,dgp_m2$res_nsel[,6],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_3_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,40), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_nsel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_nsel[,5],type="b",pch=2,lwd=3)
lines(R2s,dgp_m3$res_nsel[,6],type="b",pch=4,lwd=3)
abline(h=10,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

### Selection relevant controls

# BCCH, CV min, CV 1se

pdf("nsel_rel_1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,5), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_nsel_rel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_nsel_rel[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_nsel_rel[,3],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_rel_2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,5), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_nsel_rel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_nsel_rel[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m2$res_nsel_rel[,3],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_rel_3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,10), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_nsel_rel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_nsel_rel[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m3$res_nsel_rel[,3],type="b",pch=4,lwd=3)
abline(h=10,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

# Sensitivity to regularization choice

pdf("nsel_rel_1_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,5), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_nsel_rel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_nsel_rel[,5],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_nsel_rel[,6],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_rel_2_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,5), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_nsel_rel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_nsel_rel[,5],type="b",pch=2,lwd=3)
lines(R2s,dgp_m2$res_nsel_rel[,6],type="b",pch=4,lwd=3)
abline(h=5,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("nsel_rel_3_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,10), ylab="Number of selected controls", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_nsel_rel[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_nsel_rel[,5],type="b",pch=2,lwd=3)
lines(R2s,dgp_m3$res_nsel_rel[,6],type="b",pch=4,lwd=3)
abline(h=10,col="darkgrey",lty=2,lwd=1)
legend("right", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

############ Bias

pdf("bias_sd_1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.2), ylab="Bias/Std", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_bias[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_bias[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_bias[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_m1$res_bias[,9],type="b",pch=7,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()




### Bias/std

# BCCH, CV min, CV 1se

pdf("bias_sd_1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-0.5,2.5), ylab="Bias/Std", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,m1_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,m1_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,m1_bias_sd[,3],type="b",pch=4,lwd=3)
lines(R2s,m1_bias_sd[,9],type="b",pch=7,lwd=16)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-0.5,2.5), ylab="Bias/Std", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,m2_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,m2_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,m2_bias_sd[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-0.5,2.5), ylab="Bias/Std", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,m3_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,m3_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,m3_bias_sd[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

# Sensitivity to regularization choice

pdf("bias_sd_1_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-0.5,8), ylab="Bias/Std", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,m1_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,m1_bias_sd[,5],type="b",pch=2,lwd=3)
lines(R2s,m1_bias_sd[,6],type="b",pch=4,lwd=3)
lines(R2s,m1_bias_sd[,9],type="b",pch=7,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_2_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-0.5,8), ylab="Bias/Std", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,m2_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,m2_bias_sd[,5],type="b",pch=2,lwd=3)
lines(R2s,m2_bias_sd[,6],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_3_sens.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-0.5,8), ylab="Bias/Std", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,m3_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,m3_bias_sd[,5],type="b",pch=2,lwd=3)
lines(R2s,m3_bias_sd[,6],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topleft", c(expression(lambda[BCCH]),expression(0.5*lambda[BCCH]), expression(1.5*lambda[BCCH])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

### Coverage

pdf("cov_1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_cov[,3],type="b",pch=4,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("cov_2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m2$res_cov[,3],type="b",pch=4,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("cov_3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m3$res_cov[,3],type="b",pch=4,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

### Comparison to OLS with HCK standard error

# Coverage

pdf("cov_1_ols.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_m1$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("cov_2_ols.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m2$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_m2$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("cov_3_ols.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m3$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_m3$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

# Length

pdf("length_1_ols.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.3), ylab="Average length CI", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m1$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_m1$res_leng[,4],type="b",pch=8,lwd=3)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("length_2_ols.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.3), ylab="Average length CI", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m2$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_m2$res_leng[,4],type="b",pch=8,lwd=3)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("length_3_ols.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.3), ylab="Average length CI", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_m3$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_m3$res_leng[,4],type="b",pch=8,lwd=3)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

######################################################################################################
# Figures appendix
######################################################################################################

### Bias/std

pdf("bias_sd_app_a1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s_bern),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="DGP A1", type="n")
lines(R2s_bern,a1_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s_bern,a1_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s_bern,a1_bias_sd[,3],type="b",pch=4,lwd=3)
lines(R2s_bern,a1_bias_sd[,4],type="b",pch=8,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_app_a2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="DGP A2", type="n")
lines(R2s,a2_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,a2_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,a2_bias_sd[,3],type="b",pch=4,lwd=3)
lines(R2s,a2_bias_sd[,4],type="b",pch=8,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_app_a3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="DGP A3", type="n")
lines(R2s,a3_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,a3_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,a3_bias_sd[,3],type="b",pch=4,lwd=3)
lines(R2s,a3_bias_sd[,4],type="b",pch=8,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_app_a4.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="DGP A4", type="n")
lines(R2s,a4_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,a4_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,a4_bias_sd[,3],type="b",pch=4,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_app_a5.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="DGP A5", type="n")
lines(R2s,a5_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,a5_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,a5_bias_sd[,3],type="b",pch=4,lwd=3)
lines(R2s,a5_bias_sd[,4],type="b",pch=8,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("bias_sd_app_a6.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="DGP A6", type="n")
lines(R2s,a6_bias_sd[,1],type="b",pch=1,lwd=3)
lines(R2s,a6_bias_sd[,2],type="b",pch=2,lwd=3)
lines(R2s,a6_bias_sd[,3],type="b",pch=4,lwd=3)
lines(R2s,a6_bias_sd[,4],type="b",pch=8,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

### Coverage

pdf("cov_app_a1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s_bern),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="DGP A1", type="n")
lines(R2s_bern,dgp_a1$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s_bern,dgp_a1$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s_bern,dgp_a1$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s_bern,dgp_a1$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("cov_app_a2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="DGP A2", type="n")
lines(R2s,dgp_a2$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a2$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a2$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a2$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("cov_app_a3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="DGP A3", type="n")
lines(R2s,dgp_a3$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a3$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a3$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a3$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("cov_app_a4.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="DGP A4", type="n")
lines(R2s,dgp_a4$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a4$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a4$res_cov[,3],type="b",pch=4,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("cov_app_a5.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="DGP A5", type="n")
lines(R2s,dgp_a5$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a5$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a5$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a5$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("cov_app_a6.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0.1,1), ylab="Coverage", xlab=expression(R^{2}), main="DGP A6", type="n")
lines(R2s,dgp_a6$res_cov[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a6$res_cov[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a6$res_cov[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a6$res_cov[,4],type="b",pch=8,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

### Length

pdf("length_app_a1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s_bern),c(0,0.4), ylab="Average length CI", xlab=expression(R^{2}), main="DGP A1", type="n")
lines(R2s_bern,dgp_a1$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s_bern,dgp_a1$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s_bern,dgp_a1$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s_bern,dgp_a1$res_leng[,4],type="b",pch=8,lwd=3)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("length_app_a2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.4), ylab="Average length CI", xlab=expression(R^{2}), main="DGP A2", type="n")
lines(R2s,dgp_a2$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a2$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a2$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a2$res_leng[,4],type="b",pch=8,lwd=3)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("length_app_a3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.4), ylab="Average length CI", xlab=expression(R^{2}), main="DGP A3", type="n")
lines(R2s,dgp_a3$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a3$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a3$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a3$res_leng[,4],type="b",pch=8,lwd=3)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("length_app_a4.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.4), ylab="Average length CI", xlab=expression(R^{2}), main="DGP A4", type="n")
lines(R2s,dgp_a4$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a4$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a4$res_leng[,3],type="b",pch=4,lwd=3)
legend("bottomright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2,4), lwd=c(3,3,3), bty = "n")
graphics.off()

pdf("length_app_a5.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.4), ylab="Average length CI", xlab=expression(R^{2}), main="DGP A5", type="n")
lines(R2s,dgp_a5$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a5$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a5$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a5$res_leng[,4],type="b",pch=8,lwd=3)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()

pdf("length_app_a6.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,0.4), ylab="Average length CI", xlab=expression(R^{2}), main="DGP A6", type="n")
lines(R2s,dgp_a6$res_leng[,1],type="b",pch=1,lwd=3)
lines(R2s,dgp_a6$res_leng[,2],type="b",pch=2,lwd=3)
lines(R2s,dgp_a6$res_leng[,3],type="b",pch=4,lwd=3)
lines(R2s,dgp_a6$res_leng[,4],type="b",pch=8,lwd=3)
legend("topright", c(expression(lambda[BCCH]),expression(lambda[min]), expression(lambda["1se"]), "OLS (HCK)"), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n")
graphics.off()


######################################################################################################
# Debiased Lasso
######################################################################################################

### Bias/std

pdf("db_bias_sd_1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,m1_bias_sd[,7],type="b",pch=1,lwd=3)
lines(R2s,m1_bias_sd[,8],type="b",pch=2,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2), lwd=c(3,3), bty = "n")
graphics.off()

pdf("db_bias_sd_2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,m2_bias_sd[,7],type="b",pch=1,lwd=3)
lines(R2s,m2_bias_sd[,8],type="b",pch=2,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2), lwd=c(3,3), bty = "n")
graphics.off()

pdf("db_bias_sd_3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(-1,3), ylab="Bias/Std", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,m3_bias_sd[,7],type="b",pch=1,lwd=3)
lines(R2s,m3_bias_sd[,8],type="b",pch=2,lwd=3)
abline(h=0,col="darkgrey",lty=2,lwd=1)
legend("bottomright", c(expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2), lwd=c(3,3), bty = "n")
graphics.off()


### Coverage

pdf("db_cov_1.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,1), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=5", type="n")
lines(R2s,dgp_m1$res_cov[,7],type="b",pch=1,lwd=3)
lines(R2s,dgp_m1$res_cov[,8],type="b",pch=2,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomleft", c(expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2), lwd=c(3,3), bty = "n")
graphics.off()

pdf("db_cov_2.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,1), ylab="Coverage", xlab=expression(R^{2}), main="n=1000, p=200, k=5", type="n")
lines(R2s,dgp_m2$res_cov[,7],type="b",pch=1,lwd=3)
lines(R2s,dgp_m2$res_cov[,8],type="b",pch=2,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomleft", c(expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2), lwd=c(3,3), bty = "n")
graphics.off()

pdf("db_cov_3.pdf",pointsize=18,width=6.0,height=6.0)
plot(range(R2s),c(0,1), ylab="Coverage", xlab=expression(R^{2}), main="n=500, p=200, k=10", type="n")
lines(R2s,dgp_m3$res_cov[,7],type="b",pch=1,lwd=3)
lines(R2s,dgp_m3$res_cov[,8],type="b",pch=2,lwd=3)
abline(h=0.9,col="darkgrey",lty=2,lwd=1)
legend("bottomleft", c(expression(lambda[min]), expression(lambda["1se"])), pch=c(1,2), lwd=c(3,3), bty = "n")
graphics.off()


