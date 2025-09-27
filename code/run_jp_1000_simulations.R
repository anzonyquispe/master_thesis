#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1) Helper: the DGP and estimation functions must already be in your workspace:
#    - generate_data(model, seed)
#    - oracle_dr(dat)         → returns named c(est=..., se=...)
#    - dml_dr(dat, est, p=NULL, K=5)
#    - adml1_dr(dat, p, K=5)
#    - adml2_dr(dat, p, K=5)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
rm(list = ls())
setwd('C:/Users/eunic/Dropbox/data_eco/Udesa/mater_thesis/code')
B <- 5
set.seed(123)  # for reproducibility of the simulation seeds
source("util_functions_jp.R")
# 2) Predefine all estimator labels in the order we'll store them:
est_labels <- c(
  "theta0",
  "oracle",
  paste0("dml.ols.K", c(5,10)),
  paste0("dml.nw.K",  c(5,10)),
  paste0("dml.poly2.K", c(5,10)),
  paste0("dml.poly5.K", c(5,10)),
  paste0("dml.poly10.K",c(5,10)),
  paste0("adml1.poly2.K", c(5,10)),
  paste0("adml1.poly5.K", c(5,10)),
  paste0("adml1.poly10.K",c(5,10)),
  paste0("adml2.poly2.K", c(5,10)),
  paste0("adml2.poly5.K", c(5,10)),
  paste0("adml2.poly10.K",c(5,10))
)
P <- length(est_labels)

# 3) Storage for estimates and ses:
est_mat <- matrix(NA_real_, nrow=B, ncol=P,
                  dimnames = list(NULL, est_labels))
se_mat  <- matrix(NA_real_, nrow=B, ncol=P,
                  dimnames = list(NULL, est_labels))
theta0_vec <- numeric(B)

# 4) Simulation loop
for (b in seq_len(B)) {
  ## 4.1 draw DGP
  dat <- generate_data(model="probit", seed=b, n=400)
  theta0 <- mean(dat$gx)
  theta0_vec[b] <- theta0
  
  ## 4.2 oracle
  or <- oracle_dr(dat)
  
  ## 4.3 collect
  est_mat[b, "theta0"] <- theta0
  se_mat[b,  "theta0"] <- NA
  
  est_mat[b, "oracle"] <- or["est"]
  se_mat[b,  "oracle"] <- or["se"]
  
  ## 4.4 DML family
  i <- 3
  for (method in c("ols","nw")) {
    for (K in c(5,10)) {
      res <- dml_dr(dat, est=method, K=K)
      est_mat[b, i] <- res["est"]
      se_mat[b,  i] <- res["se"]
      i <- i+1
    }
  }
  ## 4.5 dml.poly family
  for (p in c(2,5,10)) {
    for (K in c(5,10)) {
      res <- dml_dr(dat, est="poly", p=p, K=K)
      est_mat[b, i] <- res["est"]
      se_mat[b,  i] <- res["se"]
      i <- i+1
    }
  }
  ## 4.6 ADML1
  for (p in c(2,5,10)) {
    for (K in c(5,10)) {
      res <- adml1_dr(dat, p=p, K=K)
      est_mat[b, i] <- res["est"]
      se_mat[b,  i] <- res["se"]
      i <- i+1
    }
  }
  ## 4.7 ADML2
  for (p in c(2,5,10)) {
    for (K in c(5,10)) {
      res <- adml2_dr(dat, p=p, K=K)
      est_mat[b, i] <- res["est"]
      se_mat[b,  i] <- res["se"]
      i <- i+1
    }
  }
}

# 5) Build summary table
summary_tbl <- data.frame(
  estimator   = est_labels,
  mean.theta0 = mean(theta0_vec),
  mean.est    = colMeans(est_mat, na.rm=TRUE),
  mean.bias   = colMeans( theta0_vec - est_mat , na.rm=TRUE),
  mean.se     = colMeans(se_mat, na.rm=TRUE)
)

# 6) Print nicely
print(summary_tbl, digits=3)
