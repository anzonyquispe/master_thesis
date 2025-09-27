#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Parallelized simulation for B replications
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

rm(list = ls())
setwd('/groups/sgulzar/sa_fires/proj_alterantive/adml_estimations/code/second_version/')

# 0) load your util functions
source("util_function_X_matrix.R")

# 1) setup parallel backend
library(foreach)
library(doParallel)

ncores <- parallel::detectCores() - 1
cl     <- makeCluster(ncores)
registerDoParallel(cl)

# 2) simulation settings
B <- 1000
est_labels <- c(
  "theta0",
  "oracle",
  paste0("dml.ols.K",   c(5,10)),
  paste0("dml.nw.K",    c(5,10)),
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

# 3) run the replications in parallel
res_list <- foreach(b = 1:B,
                    .packages = c("np"),             # add any pkgs your util functions need
                    .export   = c("generate_data",
                                  "oracle_dr","dml_dr",
                                  "adml1_dr","adml2_dr")
) %dopar% {
  
  tryCatch({
      
    # 3.1 draw DGP
    # dat    <- generate_data(model = "probit", seed = b, n = 400)
    dat <- generate_data_toeplitz(model = "probit", seed = b, n = 400)
    theta0 <- mean(dat$gx)
    
    inverse_true_alpha <- dat$r
    
    # 3.2 oracle
    or <- oracle_dr(dat)
    
    # 3.3 containers
    est <- numeric(P)
    se  <- numeric(P)
    alphaMSE  <- list()
    
    # 3.4 fill in
    est[1] <- theta0;    se[1] <- NA; alphaMSE[[1]] <- (inverse_true_alpha - inverse_true_alpha)
    est[2] <- or["est"]; se[2] <- or["se"]; alphaMSE[[2]] <- (inverse_true_alpha - inverse_true_alpha)
    
    
    idx <- 3
    # DML: ols & nw
    for(m in c("ols","nw")) for(K in c(5,10)) {
      tmp     <- dml_dr(dat, est = m, K = K)
      est[idx] <- tmp$est
      se[idx]  <- tmp$se
      alphaMSE[[idx]] <- (inverse_true_alpha - tmp$rest)
      idx      <- idx + 1
    }
    # DML.poly
    for(p in c(2,5,10)) for(K in c(5,10)) {
      tmp     <- dml_dr(dat, est = "poly", p = p, K = K)
      est[idx] <- tmp$est
      se[idx]  <- tmp$se
      alphaMSE[[idx]] <- (inverse_true_alpha - tmp$rest)
      idx      <- idx + 1
    }
    # ADML1
    for(p in c(2,5,10)) for(K in c(5,10)) {
      tmp     <- adml1_dr(dat, p = p, K = K)
      est[idx] <- tmp$est
      se[idx]  <- tmp$se
      alphaMSE[[idx]] <- (inverse_true_alpha - (1/tmp$alpha))
      idx      <- idx + 1
    }
    # ADML2
    for(p in c(2,5,10)) for(K in c(5,10)) {
      tmp     <- adml2_dr(dat, p = p, K = K)
      est[idx] <- tmp$est
      se[idx]  <- tmp$se
      alphaMSE[[idx]] <- (inverse_true_alpha - (1/tmp$alpha))
      idx      <- idx + 1
    }
    
    list(est = est, se = se, theta0 = theta0, alphaMSE = alphaMSE)
  }, error = function(e) {
    # If any error occurs, return NULL for this iteration:
    NULL
  })
}

# 4) stop the cluster
stopCluster(cl)

# Optionally, remove NULL entries (iterations that failed):
res_list <- Filter(Negate(is.null), res_list)

# 5) combine results into matrices
est_mat    <- t(sapply(res_list, `[[`, "est"))
se_mat     <- t(sapply(res_list, `[[`, "se"))
theta0_vec <- sapply(res_list, `[[`, "theta0")
nspecs = dim(est_mat)[2]
alphamse <- numeric(nspecs)
for (val in 1:nspecs){
  alpha.vals <- numeric(0)
  for (bsp in 1:B){
    alpha.b <- res_list[[bsp]]$alphaMSE[[val]]    
    alpha.vals <- c(alpha.vals, alpha.b)
  }
  alphamse[val] <- mean((alpha.b)**2)
}


# 6) build summary table
summary_tbl <- data.frame(
  estimator   = est_labels,
  mean.theta0 = mean(theta0_vec),
  mean.est    = colMeans(est_mat,   na.rm = TRUE),
  mean.bias   = colMeans(theta0_vec - est_mat, na.rm = TRUE),
  mean.se     = colMeans(se_mat,    na.rm = TRUE),
  alpha.mse   = alphamse
)

# 7) print
estimator_labels <- c(
  "θ₀",
  "Oracle",
  "DML (OLS, K=5)",
  "DML (OLS, K=10)",
  "DML (NW, K=5)",
  "DML (NW, K=10)",
  "DML (Poly2, K=5)",
  "DML (Poly2, K=10)",
  "DML (Poly5, K=5)",
  "DML (Poly5, K=10)",
  "DML (Poly10, K=5)",
  "DML (Poly10, K=10)",
  "ADML1 (Poly2, K=5)",
  "ADML1 (Poly2, K=10)",
  "ADML1 (Poly5, K=5)",
  "ADML1 (Poly5, K=10)",
  "ADML1 (Poly10, K=5)",
  "ADML1 (Poly10, K=10)",
  "ADML2 (Poly2, K=5)",
  "ADML2 (Poly2, K=10)",
  "ADML2 (Poly5, K=5)",
  "ADML2 (Poly5, K=10)",
  "ADML2 (Poly10, K=5)",
  "ADML2 (Poly10, K=10)"
)

summary_tbl$estimator <- estimator_labels
library(xtable)

# Format the table
latex_table <- xtable(
  summary_tbl,
  caption = "Simulation Results (Probit): Bias and Alpha MSE",
  label = "tab:sim_results",
  align = c("l", "l", rep("c", ncol(summary_tbl) - 1))
)

# Set digits for numeric columns (optional)
digits(latex_table) <- c(0, 6, 6, 6, 6, 6, 6)

# Save to .tex file
print(
  latex_table,
  file = "/groups/sgulzar/sa_fires/proj_alterantive/adml_estimations/output/second_version/simulation_results_probit_table.tex",
  include.rownames = FALSE,
  booktabs = TRUE,
  caption.placement = "top"
)
