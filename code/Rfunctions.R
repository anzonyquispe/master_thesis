######################################################################################################
# Generic functions
# Authors: Kaspar Wuthrich and Ying Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
######################################################################################################

### Generate data with pre-specified R2

gen_design <- function(n,p,k,RsqY,RsqD,alpha0,dgp="homosc", seedval = 1){
  
  set.seed(seedval)
  tbeta <- tgamma <- rep(0,p)
  tbeta[1:k] <- tgamma[1:k] <- 1
  
  # Homoscedastic dgp; allows alpha0!=0
  # if(dgp=="homosc"){
  #   
  #   cD <- sqrt(RsqD/(k-RsqD*k))
  #   
  #   if (alpha0==0){
  #     cY <- sqrt(RsqY/(k-RsqY*k))
  #   }
  #   if (alpha0!=0){
  #     a   <- k*(RsqY-1)
  #     b   <- 2*alpha0*cD*k*(RsqY-1)
  #     c   <- alpha0^2*cD^2*k*(RsqY-1)+(alpha0^2+1)*RsqY
  #     cY  <- max(Re(polyroot(c(c,b,a))))
  #   }
  # 
  #   gamma <-  cD*tgamma
  #   beta  <-  cY*tbeta
  #   
  #   X <- matrix(rnorm(n*p),n,p)
  #   d <- probit(X%*%gamma) + rnorm(n)
  #   y <- alpha0*d + X%*%beta + rnorm(n)
  #   
  # }
  # Homoscedastic dgp; allows alpha0!=0
  if(dgp=="homosc"){
    
    cD <- sqrt(RsqD/(k-RsqD*k))
    
    if (alpha0==0){
      cY <- sqrt(RsqY/(k-RsqY*k))
    }
    if (alpha0!=0){
      a   <- k*(RsqY-1)
      b   <- 2*alpha0*cD*k*(RsqY-1)
      c   <- alpha0^2*cD^2*k*(RsqY-1)+(alpha0^2+1)*RsqY
      cY  <- max(Re(polyroot(c(c,b,a))))
    }
    
    gamma <-  cD*tgamma
    beta  <-  cY*tbeta
    
    X <- matrix(rnorm(n*p),n,p)
    X <- matrix(rnorm(n*p), n, p)
    
    eta <- as.numeric(X %*% gamma)
    
    # choose target prevalence:
    pi_target <- 0.5  # change to 0.3, 0.6, etc.
    
    # find delta so that mean(plogis(eta + delta)) ~= pi_target
    f <- function(delta) mean(plogis(eta + delta)) - pi_target
    delta <- uniroot(f, c(-15, 15))$root
    
    p <- plogis(eta + delta)
    d <- rbinom(n, size = 1, prob = p)
    
    # outcome
    cY   <- sqrt(RsqY/(k - RsqY*k))
    beta <- cY * tbeta
    y <- alpha0 * d + as.numeric(X %*% beta) + rnorm(n)
  }
  
  # Heteroscedastic dgp; designed for alpha0=0
  if(dgp=="heterosc"){
    
    cD <- sqrt(RsqD/(k-RsqD*k))
    cY <- sqrt(RsqY/(k-RsqY*k))
    
    gamma <-  cD*tgamma
    beta  <-  cY*tbeta
    
    X     <- matrix(rnorm(n*p),n,p)
    sig_d <- sqrt((1+X%*%gamma)^2/mean((1+X%*%gamma)^2))
    sig_y <- sqrt((1+X%*%beta)^2/mean((1+X%*%beta)^2))
    d     <- X%*%gamma + sig_d*rnorm(n)
    y     <- X%*%beta + sig_y*rnorm(n)
    
  }
  
  # Bernoulli regressors; designed for alpha0=0
  if(dgp=="bern"){
    
    cD <- sqrt(RsqD/(k*0.25*(1-RsqD)))
    cY <- sqrt(RsqY/(k*0.25*(1-RsqY)))
    
    gamma <-  cD*tgamma
    beta  <-  cY*tbeta
    
    X <- matrix(rbinom(n*p,1,0.5),n,p)
    d <- X%*%gamma + rnorm(n)
    y <- X%*%beta + rnorm(n)
    
  }
  
  # t-distributed errors; designed for alpha0=0
  if(dgp=="tdistr"){
    
    cD <- sqrt(RsqD/(k-RsqD*k))
    cY <- sqrt(RsqY/(k-RsqY*k))
    
    gamma <-  cD*tgamma
    beta  <-  cY*tbeta
    
    X <- matrix(rnorm(n*p),n,p)
    d <- X%*%gamma + rt(n,df=5)/sqrt(5/3)
    y <- X%*%beta + rt(n,df=5)/sqrt(5/3)
    
  }
  
  return(list(X=X,d=d,y=y))
  
}

### Post double Lasso with CV regularization choice

pdl_cv <- function(y,d,X,cv){
  
  cvfit1  <- cv.glmnet(X,d,nfolds=5)
  cvfit2  <- cv.glmnet(X,y,nfolds=5)
  
  if (cv=="min"){
    gamma   <- as.matrix(coef(cvfit1, s = "lambda.min")[-1])
    pi      <- as.matrix(coef(cvfit2, s = "lambda.min")[-1])
  }
  if (cv=="1se"){
    gamma   <- as.matrix(coef(cvfit1, s = "lambda.1se")[-1])
    pi      <- as.matrix(coef(cvfit2, s = "lambda.1se")[-1])
  }
  
  ind <- (abs(gamma) + abs(pi) > 0)
  
  if (sum(ind)==0) {
    obj <- lm(y ~ d)
  }
  if (sum(ind)>0){
    Xsel  <- X[,ind]
    obj   <- lm(y ~ d + Xsel)
  } 
  
  vcov      <- vcovHC(obj, type = "HC1")
  robust_se <- sqrt(diag(vcov))
  
  alpha_hat <- obj$coefficients[2]
  se_hat    <- robust_se[2]
  
  return(list(alpha_hat=alpha_hat,se_hat=se_hat,sel_index=ind))
  
}

### lasso with DML and cross fitting
# -----------------------------------------------
# Cross-Fitted DML (Partialling-Out) with glmnet
# -----------------------------------------------
# y : numeric vector (outcome)
# d : numeric vector (treatment)
# X : matrix or data.frame of controls (will be coerced to matrix)
# K : number of folds for cross-fitting (default 5)
# R : number of repeated cross-fitting runs (>= 5 recommended)
# cv_rule : "min" or "1se" lambda choice from cv.glmnet
# seed : optional integer; if given, ensures reproducibility
#
# Returns:
#   $alpha_hat       : aggregated estimate across repeats (median)
#   $se_hat          : SE from the run whose alpha is closest to the median
#   $ci95            : 95% CI using that SE
#   $per_run         : data.frame with alpha_r, se_r, J_r, n per repeat
#   $last_run        : list with residuals and nuisances for debugging
#
# Requires: glmnet
# install.packages("glmnet")

# Requires: glmnet
# install.packages("glmnet")

dml_lasso_cf_once <- function(y, d, X, K = 5, cv_rule = c("min", "1se"), seed = NULL) {
  cv_rule <- match.arg(cv_rule)
  if (!is.null(seed)) set.seed(seed)
  
  y <- as.numeric(y); d <- as.numeric(d); X <- as.matrix(X)
  n <- length(y)
  p <- ncol(X)
  xnames <- colnames(X)
  if (is.null(xnames)) xnames <- paste0("X", seq_len(p))
  
  # build folds
  folds <- sample(rep(1:K, length.out = n))
  mhat <- ghat <- rep(NA_real_, n)
  
  # store selected indices per fold
  sel_by_fold <- vector("list", K)
  
  for (k in 1:K) {
    tr <- folds != k; te <- folds == k
    
    # ---- y | X ----
    cv_y <- glmnet::cv.glmnet(
      x = X[tr, , drop = FALSE], y = y[tr],
      nfolds = 5, family = "gaussian", standardize = TRUE
    )
    s_y <- if (cv_rule == "min") cv_y$lambda.min else cv_y$lambda.1se
    # predict OOF
    mhat[te] <- as.numeric(predict(cv_y, newx = X[te, , drop = FALSE], s = s_y))
    # selected (nonzero) controls at chosen lambda (exclude intercept)
    beta_y <- as.vector(stats::coef(cv_y, s = s_y))[-1]
    idx_y  <- which(beta_y != 0L)
    names_y <- xnames[idx_y]
    
    # ---- d | X ----
    cv_d <- glmnet::cv.glmnet(
      x = X[tr, , drop = FALSE], y = d[tr],
      nfolds = 5, family = "gaussian", standardize = TRUE
    )
    s_d <- if (cv_rule == "min") cv_d$lambda.min else cv_d$lambda.1se
    # predict OOF
    ghat[te] <- as.numeric(predict(cv_d, newx = X[te, , drop = FALSE], s = s_d))
    # selected controls for d-model
    beta_d <- as.vector(stats::coef(cv_d, s = s_d))[-1]
    idx_d  <- which(beta_d != 0L)
    names_d <- xnames[idx_d]
    
    # store selections for this fold
    sel_by_fold[[k]] <- list(
      fold = k,
      idx_y = idx_y, names_y = names_y,
      idx_d = idx_d, names_d = names_d,
      lambda_y = s_y, lambda_d = s_d
    )
  }
  
  # residualize
  u <- y - mhat
  v <- d - ghat
  
  denom <- mean(v^2)
  if (!is.finite(denom) || denom <= .Machine$double.eps)
    stop("Variance of (D - g(X)) is ~0; DML not identified with current nuisances.")
  
  # DML estimator
  alpha <- mean(v * u) / denom
  
  # IF-based SE
  psi <- v * (u - alpha * v)
  J   <- -mean(v^2)
  se  <- sqrt(mean(psi^2) / (J^2) / n)
  
  # also return union across folds if useful
  union_idx_y <- sort(unique(unlist(lapply(sel_by_fold, `[[`, "idx_y"))))
  union_idx_d <- sort(unique(unlist(lapply(sel_by_fold, `[[`, "idx_d"))))
  
  list(
    alpha_hat = alpha,
    se_hat = se,
    selected_by_fold = sel_by_fold,
    union_selected = list(
      idx_y = union_idx_y, names_y = xnames[union_idx_y],
      idx_d = union_idx_d, names_d = xnames[union_idx_d]
    )
  )
}

reg_hck <- function(y, d, X) {
  
  n  <- length(y)
  Xc <- cbind(1, X)
  dXc <- cbind(1, d, X)
  
  # First stage: d ~ X
  XXinv <- chol2inv(chol(crossprod(Xc)))
  b_d   <- XXinv %*% crossprod(Xc, d)
  res_d <- drop(d - Xc %*% b_d)
  
  # Second stage: y ~ d + X
  dXXinv <- chol2inv(chol(crossprod(dXc)))
  b_y    <- dXXinv %*% crossprod(dXc, y)
  res_y  <- drop(y - dXc %*% b_y)
  
  # Projection matrix P = Xc (Xc'Xc)^{-1} Xc'
  # Build as sparse
  P <- Xc %*% XXinv %*% t(Xc)
  P <- Matrix(P, sparse = TRUE)
  
  # M = I - P
  M <- Diagonal(n) - P
  
  # Elementwise square of M (still sparse!)
  MM <- M * M
  
  # Solve (MM) z = res_y^2 efficiently
  rhs <- res_y^2
  z   <- solve(MM, rhs, sparse = TRUE)
  
  # Now compute Sig
  Sig <- (1/n) * crossprod(res_d^2, z)
  
  Gam <- mean(res_d^2)
  alpha_hat <- b_y[2]
  se_hat    <- sqrt((1 / Gam^2) * (Sig / n))
  return(c(alpha_hat,se_hat))
}






### OLS with HCK standard errors (Cattaneo et al. (2018);
# Code is partly based on the replication material available here: https://github.com/mdcattaneo/replication-CJN_2018_JASA

# reg_hck <- function(y,d,X){
#   
#   n <- length(y)
#   
#   Xc    <- cbind(1,X)
#   dXc   <- cbind(1,d,X)
#   
#   b_d   <- chol2inv(chol(crossprod(Xc)))%*%crossprod(Xc,d)
#   res_d <- d - Xc %*% b_d
#   b_y   <- chol2inv(chol(crossprod(dXc)))%*%crossprod(dXc,y)
#   res_y <- y - dXc %*% b_y
#   
#   M     <- diag(rep(n,1))-Xc%*%chol2inv(chol(crossprod(Xc)))%*%t(Xc)
#   Sig   <- (1/n)*t(res_d^2) %*% (chol2inv(chol(M*M)))%*%res_y^2
#   Gam   <- mean(res_d^2)
#   
#   alpha_hat <- b_y[2]
#   se_hat    <- sqrt((1/Gam^2)*(Sig/n))
#   
#   return(c(alpha_hat,se_hat))
#   
# }

### OLS with robust standard errors

reg_robust <- function(y,d,X,constonly=0){
  
  if (constonly==0) {ols <- lm(y ~ d + X)}
  if (constonly==1) {ols <- lm(y ~ d)}
  
  alpha_hat   <- ols$coef[2]
  vcov        <- vcovHC(ols, type = "HC1")
  robust_se   <- sqrt(diag(vcov))
  se_hat      <- robust_se[2]
  
  return(list(alpha_hat=alpha_hat,se_hat=se_hat))
}

### Debiased Lasso (van der Geer et al. (2014)); code designed for homoscedastic errors
# See Caner & Kock (2018) for an extension and code for heteroscedastic errors

db_cv <- function(y,d,X,cv){
  
  n <- length(y)
  
  y <- y - mean(y)
  d <- d - mean(d)
  X <- scale(X, scale = FALSE)
  
  cvfit1  <- cv.glmnet(X,d,nfolds=5,intercept=FALSE)
  cvfit2  <- cv.glmnet(cbind(d,X),y,nfolds=5,intercept=FALSE)
  
  if (cv=="min"){
    gamma   <- as.matrix(coef(cvfit1, s = "lambda.min")[-1])
    beta    <- as.matrix(coef(cvfit2, s = "lambda.min")[-1])
  }
  if (cv=="1se"){
    gamma   <- as.matrix(coef(cvfit1, s = "lambda.1se")[-1])
    beta    <- as.matrix(coef(cvfit2, s = "lambda.1se")[-1])
  }
  
  tau2      <- (t(d - X %*% gamma) %*% d)/n
  Theta     <- c(1/tau2)*c(1,-gamma)
  bias      <- (Theta %*% t(cbind(d,X)) %*% (y-cbind(d,X) %*% beta))/n
  alpha_hat <- beta[1] + bias
  u_hat     <- c(y - cbind(d,X)%*%beta)
  se_hat    <- sqrt(mean(u_hat^2))*sqrt((Theta %*% ((t(cbind(d,X)) %*%  cbind(d,X))/n) %*% Theta)/n)

  return(list(alpha_hat=alpha_hat,se_hat=se_hat))
  
}

