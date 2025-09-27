# R translation of the MATLAB simulation using random forests and sample‐splitting
rm(list = ls())
# Load required packages
library(MASS)           # for mvrnorm
library(foreach)        # for parallel foreach
library(doParallel)     # for parallel backend
library(glmnet)
setwd("C:/Users/Anzony/Dropbox/data_eco/Udesa/mater_thesis")
source("code/0_riesz_code_for_simulations.R")

# Simulation parameters
set.seed(2162016)
nrep   <- 5000
n      <- 500
p      <- 20
alpha  <- 0.5

# Construct covariance for X ~ N(0, Σ) with Σ_{ij} = 0.7^{|i-j|}
Sigma <- toeplitz(0.7^(0:(p-1)))

# Preallocate result matrices
# Columns: 1 = estimate, 2 = standard error
rfsss1  <- matrix(0, nrep, 2)
rfgsss1 <- matrix(0, nrep, 2)
rfsds1  <- matrix(0, nrep, 2)
rfsss2  <- matrix(0, nrep, 2)
rfgsss2 <- matrix(0, nrep, 2)
rfsds2  <- matrix(0, nrep, 2)
rfsss   <- matrix(0, nrep, 2)
rfgsss  <- matrix(0, nrep, 2)
rfsds   <- matrix(0, nrep, 2)

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Run the simulation in parallel
res_list <- foreach(ii = 1:nrep, .combine = rbind, .packages = c("MASS","glmnet", "nnet", "hdm")) %dopar% {
  
  # 1) Draw covariates
  set.seed(ii)
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  # 2) Generate D and Y
  # d <- 1 * X[,1] + 0.25 * exp(X[,3])/(1 + exp(X[,3])) + rnorm(n, sd = 1)
  lin_pred <- 1 * X[, 1] + 0.25 * exp(X[, 3]) / (1 + exp(X[, 3])) + rnorm(n, sd = 1)
  prob_D   <- pnorm(lin_pred)          # convert to probability via logistic link
  d        <- rbinom(n, size = 1, prob = prob_D)
  
  y <- alpha * d + 1 * exp(X[,1])/(1 + exp(X[,1])) + 0.25 * X[,3] + rnorm(n, sd = 1)
  
  # 3) Sample splitting
  train_idx <- sample.int(n, floor(n/2))
  test_idx  <- setdiff(seq_len(n), train_idx)
  folds <- list(
    `1` = train_idx,
    `2` = test_idx
  )
  
  
  ## First split ##
  # (a) Predict Y on test set
  fit_y1   <- cv.glmnet(X[train_idx,], y[train_idx], alpha = 1)
  fit_y1   <- glmnet(X[train_idx,], y[train_idx], alpha = 1, lambda=fit_y1$lambda.min)
  yhat1    <- predict(fit_y1, X[test_idx,])
  ry1      <- y[test_idx] - yhat1
  
  # (b) Predict D on test set
  fit_d1   <- cv.glmnet(X[train_idx,], d[train_idx], alpha = 1)
  fit_d1   <- glmnet(X[train_idx,], d[train_idx], alpha = 1, lambda=fit_d1$lambda.min)
  dhat1    <- predict(fit_d1, X[test_idx,])
  rd1      <- d[test_idx] - dhat1
  
  # (c) Non‐Orthogonal (D) + sample‐splitting
  beta_sss1 <- sum(rd1 * y[test_idx]) / sum(rd1 * d[test_idx])
  se_sss1   <- sqrt(
    mean((y[test_idx] - d[test_idx] * beta_sss1)^2) *
      sum(rd1^2) /
      (sum(rd1 * d[test_idx])^2)
  )
  
  # (d) Non‐Orthogonal (Y) + sample‐splitting
  ghat1     <- yhat1 - dhat1 * beta_sss1
  gy1       <- y[test_idx] - ghat1
  beta_gsss1 <- sum(d[test_idx] * gy1) / sum(d[test_idx]^2)
  se_gsss1   <- sqrt(
    sum((gy1 - d[test_idx] * beta_gsss1)^2) / (length(gy1) - 1) *
      (1 / sum(d[test_idx]^2))
  )
  
  # (e) Orthogonal + sample‐splitting (direct residual regression)
  beta_sds1 <- sum(rd1 * ry1) / sum(rd1^2)
  se_sds1   <- sqrt(
    sum((ry1 - rd1 * beta_sds1)^2) / (length(ry1) - 1) *
      (1 / sum(rd1^2))
  )
  
  ## Second split ##
  fit_y2   <- cv.glmnet(X[test_idx,], y[test_idx], alpha = 1)
  fit_y2   <- glmnet(X[test_idx,], y[test_idx], alpha = 1, lambda=fit_y2$lambda.min)
  yhat2    <- predict(fit_y2, X[train_idx,])
  ry2      <- y[train_idx] - yhat2
  
  fit_d2   <- cv.glmnet(X[test_idx,], d[test_idx], alpha = 1)
  fit_d2   <- glmnet(X[test_idx,], d[test_idx], alpha = 1, lambda=fit_d2$lambda.min)
  dhat2    <- predict(fit_d2, X[train_idx,])
  rd2      <- d[train_idx] - dhat2
  
  beta_sss2 <- sum(rd2 * y[train_idx]) / sum(rd2 * d[train_idx])
  se_sss2   <- sqrt(
    mean((y[train_idx] - d[train_idx] * beta_sss2)^2) *
      sum(rd2^2) /
      (sum(rd2 * d[train_idx])^2)
  )
  
  ghat2     <- yhat2 - dhat2 * beta_sss2
  gy2       <- y[train_idx] - ghat2
  beta_gsss2 <- sum(d[train_idx] * gy2) / sum(d[train_idx]^2)
  se_gsss2   <- sqrt(
    sum((gy2 - d[train_idx] * beta_gsss2)^2) / (length(gy2) - 1) *
      (1 / sum(d[train_idx]^2))
  )
  
  beta_sds2 <- sum(rd2 * ry2) / sum(rd2^2)
  se_sds2   <- sqrt(
    sum((ry2 - rd2 * beta_sds2)^2) / (length(ry2) - 1) *
      (1 / sum(rd2^2))
  )
  
  ## Aggregate splits ##
  beta_sss  <- 0.5 * (beta_sss1 + beta_sss2)
  se_sss    <- sqrt(0.25 * se_sss1^2 + 0.25 * se_sss2^2)
  beta_gsss <- 0.5 * (beta_gsss1 + beta_gsss2)
  se_gsss   <- sqrt(0.25 * se_gsss1^2 + 0.25 * se_gsss2^2)
  beta_sds  <- 0.5 * (beta_sds1 + beta_sds2)
  se_sds    <- sqrt(0.25 * se_sds1^2 + 0.25 * se_sds2^2)
  
  res.adml.dml <- adml_res(y, d, X)
  beta_adml <- res.adml.dml[1]
  se_adml <- res.adml.dml[2]
  beta_dml <- res.adml.dml[3]
  se_dml <- res.adml.dml[4]
  
  c(beta_sss, se_sss, beta_gsss, se_gsss, beta_sds, se_sds, beta_adml, se_adml, beta_dml, se_dml)
}

# Stop parallel cluster
stopCluster(cl)

# Reshape results
res_mat    <- do.call(rbind, lapply(seq_len(nrow(res_list)), function(i) res_list[i,]))
rfsss      <- res_mat[, 1:2]
rfgsss     <- res_mat[, 3:4]
rfsds      <- res_mat[, 5:6]
rfadml      <- res_mat[, 7:8]
rfdml      <- res_mat[, 9:10]


# Coverage calculations


cov_sss  <- mean((rfsss[,1]-1.96*rfsss[,2]<=alpha) & (rfsss[,1]+1.96*rfsss[,2]>=alpha))
cov_gsss <- mean((rfgsss[,1]-1.96*rfgsss[,2]<=alpha) & (rfgsss[,1]+1.96*rfgsss[,2]>=alpha))
cov_sds  <- mean((rfsds[,1]-1.96*rfsds[,2]<=alpha) & (rfsds[,1]+1.96*rfsds[,2]>=alpha))
cov_adml <- mean((rfadml[,1]-1.96*rfadml[,2]<=alpha) & (rfadml[,1]+1.96*rfadml[,2]>=alpha))
cov_dml  <- mean((rfdml[,1]-1.96*rfdml[,2]<=alpha) & (rfdml[,1]+1.96*rfdml[,2]>=alpha))

cat("RF - Nonorthogonal (D) + sample splitting:", cov_sss, "\n")
cat("RF - Nonorthogonal (Y) + sample splitting:", cov_gsss, "\n")
cat("RF - Orthogonal + sample splitting:", cov_sds, "\n")
cat("Lasso ADML - Orthogonal + sample splitting:", cov_adml, "\n")
cat("Lasso DML - Orthogonal + sample splitting:", cov_dml, "\n")


mean((rfsss[,1]  - alpha) )^2

bias_sds  <- mean((rfsds[,1]  - alpha) )^2
bias_adml  <- mean((rfadml[,1]  - alpha))^2


var_adml <- var(rfadml[,1])
var_sds  <- var(rfsds[,1])

mse_sss  <- mean((rfsss[,1]  - alpha) ^2)
mse_gsss <- mean((rfgsss[,1] - alpha) ^2)
mse_sds  <- mean((rfsds[,1]  - alpha) ^2)
mse_adml <- mean((rfadml[,1]  - alpha)^2)
mse_dml  <- mean((rfdml[,1]  - alpha) ^2)

cat("RF - Nonorthogonal (D) + sample splitting:", mse_sss, "\n")
cat("RF - Nonorthogonal (Y) + sample splitting:", mse_gsss, "\n")
cat("RF - Orthogonal + sample splitting:", mse_sds, "\n")
cat("Lasso ADML - Orthogonal + sample splitting:", mse_adml, "\n")
cat("Lasso DML - Orthogonal + sample splitting:", mse_dml, "\n")


cat("Lasso - Orthogonal + sample splitting:", cov_sds, "\n")
cat("Lasso ADML - Orthogonal + sample splitting:", cov_adml, "\n")
cat("Lasso - Orthogonal + sample splitting:", mse_sds, "\n")
cat("Lasso ADML - Orthogonal + sample splitting:", mse_adml, "\n")



# Set up JPEG output
jpeg("output/theta_dist.jpeg", width = 1200, height = 800)

# Arrange plots in 1 row, 2 columns
par(mfrow = c(1, 1))

edges <- seq(0, 1, by = 0.015)
hist(rfadml[,1], breaks = edges, freq = FALSE,
     col = rgb(0,0,1,0.25), main = "Simulated Distributions of θ̂", xlab = "")
hist(rfsds[,1], breaks = edges, freq = FALSE,
     col = rgb(1,0,0,0.25), add = TRUE)
abline(v = alpha, lwd = 2)
legend("topright", legend = c("ADML", "DML"),
       fill = c(rgb(0,0,1,0.25), rgb(1,0,0,0.25)))

# Close the JPEG device
dev.off()





# Plotting overlaid histograms of t-statistics
se1 <- sd(rfgsss[,1])
se2 <- sd(rfsds[,1])
se3 <- sd(rfadml[,1])
se4 <- sd(rfdml[,1])
tse1 <- se1 * (rfgsss[,1] - alpha) / rfgsss[,2]
tse2 <- se2 * (rfsds[,1] - alpha) / rfsds[,2]
tse3 <- se3 * (rfadml[,1] - alpha) / rfadml[,2]
tse4 <- se4 * (rfdml[,1] - alpha) / rfdml[,2]

# Non-Orthogonal t-stat
hist(tse1,
     breaks = floor(sqrt(nrep)), freq = FALSE,
     col = rgb(0,0.5,0.75,0.5),
     main = sprintf("Non-Orthogonal, n = %d, p = %d", n, p),
     xlab = "t-statistic")
curve(dnorm(x, mean = 0, sd = se1), add = TRUE, col = "red")
legend("topright", legend = c("Simulation", "N(0,Σ)"), fill = c(rgb(0,0.5,0.75,0.5), NA),
       border = NA, lty = c(NA,1), col = c(NA, "red"))

# Orthogonal t-stat
hist(tse2,
     breaks = floor(sqrt(nrep)), freq = FALSE,
     col = rgb(0,0.5,0.75,0.5),
     main = sprintf("Orthogonal, n = %d, p = %d", n, p),
     xlab = "t-statistic")
curve(dnorm(x, mean = 0, sd = se2), add = TRUE, col = "red")
legend("topright", legend = c("Simulation", "N(0,Σ) - DML"), fill = c(rgb(0,0.5,0.75,0.5), NA),
       border = NA, lty = c(NA,1), col = c(NA, "red"))




hist(tse3,
     breaks = floor(sqrt(nrep)), freq = FALSE,
     col = rgb(0,0.5,0.75,0.5),
     main = sprintf("Orthogonal, n = %d, p = %d", n, p),
     xlab = "t-statistic")
curve(dnorm(x, mean = 0, sd = se3), add = TRUE, col = "red")
legend("topright", legend = c("Simulation", "N(0,Σ) - ADML"), fill = c(rgb(0,0.5,0.75,0.5), NA),
       border = NA, lty = c(NA,1), col = c(NA, "red"))





# Set up JPEG output
jpeg("output/orthogonal_tstat_plots.jpeg", width = 1200, height = 600)

# Arrange plots in 1 row, 2 columns
par(mfrow = c(1, 2))

# Plot 1
hist(tse2,
     breaks = floor(sqrt(nrep)), freq = FALSE,
     col = rgb(0,0.5,0.75,0.5),
     main = sprintf("Orthogonal, n = %d, p = %d", n, p),
     xlab = "t-statistic")
curve(dnorm(x, mean = 0, sd = se2), add = TRUE, col = "red")
legend("topright", legend = c("Simulation", "N(0,Σ) - DML"),
       fill = c(rgb(0,0.5,0.75,0.5), NA), border = NA,
       lty = c(NA,1), col = c(NA, "red"))

# Plot 2
hist(tse3,
     breaks = floor(sqrt(nrep)), freq = FALSE,
     col = rgb(0,0.5,0.75,0.5),
     main = sprintf("Orthogonal, n = %d, p = %d", n, p),
     xlab = "t-statistic")
curve(dnorm(x, mean = 0, sd = se3), add = TRUE, col = "red")
legend("topright", legend = c("Simulation", "N(0,Σ) - ADML"),
       fill = c(rgb(0,0.5,0.75,0.5), NA), border = NA,
       lty = c(NA,1), col = c(NA, "red"))

# Close the JPEG device
dev.off()




# Set up JPEG output
jpeg("output/theta_dist.jpeg", width = 1200, height = 600)

# Arrange plots in 1 row, 2 columns
par(mfrow = c(1, 1))

edges <- seq(0, 1, by = 0.015)
hist(rfadml[,1], breaks = edges, freq = FALSE,
     col = rgb(0,0,1,0.25), main = "Simulated Distributions of θ̂", xlab = "")
hist(rfsds[,1], breaks = edges, freq = FALSE,
     col = rgb(1,0,0,0.25), add = TRUE)
abline(v = alpha, lwd = 2)
legend("topright", legend = c("ADML", "DML"),
       fill = c(rgb(0,0,1,0.25), rgb(1,0,0,0.25)))

# Close the JPEG device
dev.off()


# Simulated distributions of θ̂
edges <- seq(0, 1, by = 0.015)
hist(rfgsss[,1], breaks = edges, freq = FALSE,
     col = rgb(0,0,1,0.25), main = "Simulated Distributions of θ̂", xlab = "")
hist(rfsds[,1], breaks = edges, freq = FALSE,
     col = rgb(1,0,0,0.25), add = TRUE)
abline(v = alpha, lwd = 2)
legend("topright", legend = c("Non-Orthogonal", "Orthogonal"),
       fill = c(rgb(0,0,1,0.25), rgb(1,0,0,0.25)))






