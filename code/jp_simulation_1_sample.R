rm(list=ls())

B       <- 1000
## 2a.  Design
dat<-generate_data(model = "probit", seed = 1)
names(dat)

############################################################
## 0.  Parameter value θ0 
############################################################

theta0 <- mean(dat$gx)

############################################################
## 1.  Oracle-DR estimator of θ0 = E[Y]
############################################################


oracle.res <- oracle_dr(dat)

############################################################
## 2.  DML estimator of θ0
############################################################
dml.res1 <- dml_dr(dat, est='ols')
dml.res2 <- dml_dr(dat, est='nw')
dml.res3 <- dml_dr(dat, est='poly', p = 2)
dml.res4 <- dml_dr(dat, est='poly', p = 5)
dml.res5 <- dml_dr(dat, est='poly', p = 10)


dml.res6 <- dml_dr(dat, est='ols', K=10)
dml.res7 <- dml_dr(dat, est='nw', K=10)
dml.res8 <- dml_dr(dat, est='poly', p = 2, K=10)
dml.res9 <- dml_dr(dat, est='poly', p = 5, K=10)
dml.res10 <- dml_dr(dat, est='poly', p = 10, K=10)

############################################################
## 3.  ADML1 estimator of θ0
############################################################

adml1.res1 <- adml1_dr(dat, p = 2)
adml1.res2 <- adml1_dr(dat, p = 5)
adml1.res3 <- adml1_dr(dat, p = 10)

adml1.res3 <- adml1_dr(dat, p = 2, K=10)
adml1.res4 <- adml1_dr(dat, p = 5, K=10)
adml1.res5 <- adml1_dr(dat, p = 10, K=10)


############################################################
## 4.  ADML2 estimator of θ0
############################################################

adml2.res1 <- adml2_dr(dat, p = 2)
adml2.res2 <- adml2_dr(dat, p = 5)
adml2.res3 <- adml2_dr(dat, p = 10)

adml2.res3 <- adml2_dr(dat, p = 2, K=10)
adml2.res4 <- adml2_dr(dat, p = 5, K=10)
adml2.res5 <- adml2_dr(dat, p = 10, K=10)
