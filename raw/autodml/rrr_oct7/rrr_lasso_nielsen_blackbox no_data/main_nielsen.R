########
# set up
########

rm(list=ls())

library("foreign")
library("ggplot2")
library("quantreg")

library("glmnet")
library("nnet")
library("randomForest")

library("dplyr")
library("tibble")
library("tidyr")
library("MGLM")
library("multiwayvcov")

setwd("~/Documents/research/rrr_lasso_nielsen_blackbox")

#######################
# clean and format data
#######################

source('data_cleaning.R') #done
source('specifications.R') #in progress

milk_indicator=1 #0: soda, 1: milk
poly_degree=4 #degree of own price, income, interactions thereof
centering=1 #0: original data, 1: indirect componenets centered + scaled

data<-get_data(df,milk_indicator,poly_degree,centering) #specific to poly_degree=4 for grep

hh=data[[1]]
t=data[[2]]
Ti=data[[3]]
Y=data[[4]]
regressors=data[[5]]
dregressors=data[[6]]

N=length(unique(hh))
sum_Ti=length(hh)
p=dim(regressors)[2]

##################
# helper functions
##################

source('primitives.R')
source('stage1.R')

#p0 used in low-dim dictionary in the stage 1 tuning procedure
p0=10

D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
D_add=0.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights

###########
# algorithm
###########

set.seed(1) # for sample splitting

alpha_estimator=1
gamma_estimator=1
bias=0 

#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso. only, since differentiable
#bias: 0 DML, 1 plug-in

# for debug of stage 2
#hh<-hh[1:1000]
#Y<-Y[1:1000]
#X<-regressors[1:1000,]
#dX<-dregressors[1:1000,]

source('stage2.R')
results<-rrr(hh,Y,regressors,dregressors,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator,bias)
printer(results)
for_tex(results)

#}