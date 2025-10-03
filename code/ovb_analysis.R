######################################################################################################
# Code for reproducing the simulation results in the paper and appendix
# Authors: Kaspar Wuthrich and Ying Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
######################################################################################################

rm(list = ls())

setwd("C:/Users/187697/Documents/GitHub/master_thesis")

######################################################################################################


R2s       <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)
# dgp_m1 <- readRDS("data/dgp_m1.rds")
# dgp_m2 <- readRDS("data/dgp_m2.rds")
# dgp_m3 <- readRDS("data/dgp_m3.rds")

dgp_m1 <- readRDS("output/dgp_m1_adml.rds")
dgp_m2 <- readRDS("output/dgp_m2_adml.rds")
dgp_m3 <- readRDS("output/dgp_m3_adml.rds")


png("output/coverage_plotn500.png", width=800, height=600, res=120)
plot(range(R2s), c(0,1),
     ylab="Coverage",
     xlab=expression(R^{2}),
     main="p=200, k=5, n = 500",
     type="n")

# Ols
lines(R2s, dgp_m1$res_cov[,4],  type="b", pch=1, lwd=3, col="blue")
# lines(R2s, dgp_m2$res_cov[,4],  type="b", pch=2, lwd=3, col="blue")
# lines(R2s, dgp_m3$res_cov[,4],  type="b", pch=3, lwd=3, col="blue")

# DoubleLasso
lines(R2s, dgp_m1$res_cov[,1],  type="b", pch=1, lwd=3, col="pink")
# lines(R2s, dgp_m2$res_cov[,1],  type="b", pch=2, lwd=3, col="pink")
# lines(R2s, dgp_m3$res_cov[,1],  type="b", pch=3, lwd=3, col="pink")

# DML
lines(R2s, dgp_m1$res_cov[,9], type="b", pch=1, lwd=3, col="green")
# lines(R2s, dgp_m2$res_cov[,9], type="b", pch=2, lwd=3, col="green")
# lines(R2s, dgp_m3$res_cov[,9], type="b", pch=3, lwd=3, col="green")

# ADML
lines(R2s, dgp_m1$res_cov[,10], type="b", pch=1, lwd=3, col="red")
# lines(R2s, dgp_m2$res_cov[,10], type="b", pch=2, lwd=3, col="red")
# lines(R2s, dgp_m3$res_cov[,10], type="b", pch=3, lwd=3, col="red")

# Legend for pch
# legend("topright",
#        legend=c("n=500", "n=1000", "n=5000"),
#        pch=c(1, 2, 3),
#        bty="n")
legend("bottomleft",
       legend=c("OLS", "DML Lasso", "Double Lasso"),
       col=c("blue", "red", "green"),
       lwd=3,
       bty="n")
dev.off()



png("output/coverage_plotn1000.png", width=800, height=600, res=120)
plot(range(R2s), c(0,1),
     ylab="Coverage",
     xlab=expression(R^{2}),
     main="p=200, k=5, n = 1000",
     type="n")

# Ols
# lines(R2s, dgp_m1$res_cov[,4],  type="b", pch=1, lwd=3, col="blue")
lines(R2s, dgp_m2$res_cov[,4],  type="b", pch=2, lwd=3, col="blue")
# lines(R2s, dgp_m3$res_cov[,4],  type="b", pch=3, lwd=3, col="blue")

# DoubleLasso
# lines(R2s, dgp_m1$res_cov[,1],  type="b", pch=1, lwd=3, col="pink")
lines(R2s, dgp_m2$res_cov[,1],  type="b", pch=2, lwd=3, col="pink")
# lines(R2s, dgp_m3$res_cov[,1],  type="b", pch=3, lwd=3, col="pink")

# DML
# lines(R2s, dgp_m1$res_cov[,9], type="b", pch=1, lwd=3, col="green")
lines(R2s, dgp_m2$res_cov[,9], type="b", pch=2, lwd=3, col="green")
# lines(R2s, dgp_m3$res_cov[,9], type="b", pch=3, lwd=3, col="green")

# ADML
# lines(R2s, dgp_m1$res_cov[,10], type="b", pch=1, lwd=3, col="red")
lines(R2s, dgp_m2$res_cov[,10], type="b", pch=2, lwd=3, col="red")
# lines(R2s, dgp_m3$res_cov[,10], type="b", pch=3, lwd=3, col="red")

# Legend for pch
# legend("topright",
#        legend=c("n=500", "n=1000", "n=5000"),
#        pch=c(1, 2, 3),
#        bty="n")
legend("bottomleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()








png("output/coverage_plotn5000.png", width=800, height=600, res=120)
plot(range(R2s), c(0.5,0.96),
     ylab="Coverage",
     xlab=expression(R^{2}),
     main="p=200, k=5, n = 5000",
     type="n")

# Ols
lines(R2s, dgp_m3$res_cov[,4],  type="b", pch=3, lwd=3, col="blue")

# DoubleLasso
lines(R2s, dgp_m3$res_cov[,1],  type="b", pch=3, lwd=3, col="pink")

# DML
lines(R2s, dgp_m3$res_cov[,9], type="b", pch=3, lwd=3, col="green")

# ADML
lines(R2s, dgp_m3$res_cov[,10], type="b", pch=3, lwd=3, col="red")

# Legend for pch
# legend("topright",
#        legend=c("n=500", "n=1000", "n=5000"),
#        pch=c(1, 2, 3),
#        bty="n")
legend("bottomleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()



################################################################################

################################################################################
################### Bias Plot ##################################################

png("output/bias_plotn500.png", width=800, height=600, res=120)
plot(range(R2s), c(0.07,0.2),
     ylab="Bias",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=500",
     type="n")


# Ols
lines(R2s, dgp_m1$res_bias[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m1$res_bias[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m1$res_bias[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m1$res_bias[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()





png("output/bias_plotn1000.png", width=800, height=600, res=120)
plot(range(R2s), c(0.045,0.1),
     ylab="Bias",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=1000",
     type="n")


# Ols
lines(R2s, dgp_m2$res_bias[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m2$res_bias[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m2$res_bias[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m2$res_bias[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()




png("output/bias_plotn5000.png", width=800, height=600, res=120)
plot(range(R2s), c(0.02,0.031),
     ylab="Bias",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=5000",
     type="n")


# Ols
lines(R2s, dgp_m3$res_bias[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m3$res_bias[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m3$res_bias[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m3$res_bias[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()

################################################################################



################################################################################
################### Bias Plot ##################################################

png("output/sd_plotn500.png", width=800, height=600, res=120)
plot(range(R2s), c(0.08,0.15),
     ylab="SE",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=500",
     type="n")


# Ols
lines(R2s, dgp_m1$res_sd[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m1$res_sd[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m1$res_sd[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m1$res_sd[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()





png("output/sd_plotn1000.png", width=800, height=600, res=120)
plot(range(R2s), c(0.06,0.09),
     ylab="SE",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=1000",
     type="n")


# Ols
lines(R2s, dgp_m2$res_sd[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m2$res_sd[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m2$res_sd[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m2$res_sd[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()




png("output/sd_plotn5000.png", width=800, height=600, res=120)
plot(range(R2s), c(0.0275,0.032),
     ylab="SE",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=1000",
     type="n")


# Ols
lines(R2s, dgp_m3$res_sd[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m3$res_sd[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m3$res_sd[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m3$res_sd[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()

################################################################################



################################################################################
################### Bias/SE Plot ##################################################


png("output/biasse_plotn500.png", width=800, height=600, res=120)
plot(range(R2s), c(0.5,2),
     ylab="Bias/SE",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=500",
     type="n")


# Ols
lines(R2s, dgp_m1$res_bias[,4]/dgp_m1$res_sd[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m1$res_bias[,1]/dgp_m1$res_sd[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m1$res_bias[,9]/dgp_m1$res_sd[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m1$res_bias[,10]/dgp_m1$res_sd[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()





png("output/biasse_plotn1000.png", width=800, height=600, res=120)
plot(range(R2s), c(0.75,1.5),
     ylab="Bias/SE",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=1000",
     type="n")


# Ols
lines(R2s, dgp_m2$res_bias[,4]/dgp_m2$res_sd[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m2$res_bias[,1]/dgp_m2$res_sd[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m2$res_bias[,9]/dgp_m2$res_sd[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m2$res_bias[,10]/dgp_m2$res_sd[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()




png("output/biasse_plotn5000.png", width=800, height=600, res=120)
plot(range(R2s), c(0.75,1),
     ylab="Bias/SE",
     xlab=expression(R^{2}),
     main="p=200, k=5, n=1000",
     type="n")


# Ols
lines(R2s, dgp_m3$res_bias[,4]/dgp_m3$res_sd[,4],  type="b", pch=1, lwd=3, col="blue")
# Double Lasso
lines(R2s, dgp_m3$res_bias[,1]/dgp_m3$res_sd[,1],  type="b", pch=1, lwd=3, col="pink")
# DML
lines(R2s, dgp_m3$res_bias[,9]/dgp_m3$res_sd[,9], type="b", pch=1, lwd=3, col="green")
# ADML
lines(R2s, dgp_m3$res_bias[,10]/dgp_m3$res_sd[,10], type="b", pch=1, lwd=3, col="red")

# Legend for pch
legend("topleft",
       legend=c("OLS", "Double Lasso", "DML Lasso", "ADML Lasso" ),
       col=c("blue", "pink", "green", "red" ),
       lwd=3,
       bty="n")
dev.off()



################################################################################
