library(np)
library("matrixcalc")

#### ADML function
adml <- function(dat) {
  gamma <- dat$gx        # true γ0(x)
  alpha <- dat$alpha           # true α0(x)
  alpha * dat$D * (dat$Y - gamma) + gamma 
}

adml2_dr <- function(dat,p=2, K = 5 ) {
  
  
  D <- dat$D
  X <- dat$X
  n <- dim(X)[1]
  Y <- dat$Y
  gx <- dat$gx
  alpha0 <- 1/dat$r
  
  
  fold_size <- n / K
  Psi_tilde.adml2 <- numeric(0)
  alpha.est <- numeric(0)
  
  X <- as.matrix(X)
  d <- ncol(X)
  BX <- matrix(1, n, 1) 
  for (j in 1:d) {
    for (deg in 1:p) {
      BX <- cbind(BX, X[, j]^deg)
    }
  }
  
  
  for (k in 1:K) {
    idx_start <- (k - 1) * fold_size + 1
    idx_end <- k * fold_size
    idx <- idx_start:idx_end
    
    # Select the test fold
    d.sel <- D[idx]
    b.sel <- BX[idx, ]
    y.sel <- Y[idx]
    
    # Alpha hat estimation using OLS
    # Compute each component
    # Subset the data
    D_sub <- D[-idx]
    bX_sub <- BX[-idx, , drop = FALSE]
    alpha0_sub <- alpha0[-idx]
    
    W <- bX_sub * D_sub  # element-wise multiply each row by corresponding D_i
    A <- t(W) %*% bX_sub  # (p+1) x (p+1) matrix
    B <- t(W) %*% alpha0_sub  # (p+1) x 1 vector
    
    # Compute gamma_hat
    alpha.coeff <- solve(A, B)  # (p+1) x 1 vector
    alpha_hat <- b.sel %*% alpha.coeff
    # Saving alpha values
    alpha.est <- c(alpha.est, alpha_hat)
    
    # Appending data frame using the main columns
    dat <- data.frame(X = b.sel, D = d.sel, Y = y.sel, alpha = alpha_hat, gx = gx[idx])
    Psi_tilde.l.adml2 <- adml(dat)
    
    # Store the result
    Psi_tilde.adml2 <- c(Psi_tilde.adml2, Psi_tilde.l.adml2)
  }
  theta.adml2 <- mean(Psi_tilde.adml2)
  devs <- Psi_tilde.adml2-theta.adml2
  se.adml2 <- sqrt(mean((devs)^2)/n)
  return(list("est" =theta.adml2, 
              "se" = se.adml2,
              "alpha" = alpha.est))
}



adml1_dr <- function(dat,p=2, K = 5 ) {
  
  D <- dat$D
  X <- dat$X
  Y <- dat$Y
  gx <- dat$gx
  
  
  # Define the indices for the k-th fold
  n <- dim(X)[1]
  fold_size <- n / K
  Psi_tilde.adml1 <- numeric(0)
  alpha.est <- numeric(0)
  
  X <- as.matrix(X)
  d <- ncol(X)
  BX <- matrix(1, n, 1) 
  for (j in 1:d) {
    for (deg in 1:p) {
      BX <- cbind(BX, X[, j]^deg)
    }
  }
  
  for (k in 1:K) {
    idx_start <- (k - 1) * fold_size + 1
    idx_end <- k * fold_size
    idx <- idx_start:idx_end
    
    # Select the test fold
    d.sel <- D[idx]
    b.sel <- BX[idx, ]
    y.sel <- Y[idx]
    
    # Subset
    D_sub <- D[-idx]
    bX_sub <- BX[-idx,]
    
    # Weighted matrix (each row of bX_sub is scaled by D_i)
    W <- bX_sub * D_sub
    
    # Compute matrix A and vector B
    A <- t(W) %*% bX_sub  # (p+1) x (p+1)
    B <- colSums(bX_sub)  # (p+1) x 1
    
    # Compute gamma_hat
    alpha.coeff <- solve(A, B)  # (p+1) x 1 vector
    alpha_hat <- b.sel %*% alpha.coeff
    
    # Saving alpha values
    alpha.est <- c(alpha.est, alpha_hat)
    # Appending data frame using the main columns
    dat <- data.frame(X = b.sel, D = d.sel, Y = y.sel, alpha = alpha_hat, gx = gx[idx])
    Psi_tilde.l.adml1 <- adml(dat)
    
    # Store the result
    Psi_tilde.adml1 <- c(Psi_tilde.adml1, Psi_tilde.l.adml1)
  }
  theta.adml1 <- mean(Psi_tilde.adml1)
  devs <- Psi_tilde.adml1-theta.adml1
  se.adml1 <- sqrt(mean((devs)^2)/n)
  return(list("est" =theta.adml1, 
              "se" = se.adml1, 
              "alpha" = alpha.est))
}

dml_dr <- function(dat,est='ols',p=2, K = 5) {
  
  
  D <- as.matrix(dat$D)
  n <- dim(D)[1]
  X <- as.matrix(dat$X)
  Y <- dat$Y
  gx <- dat$gx
  
  fold_size <- n / K
  Psi_tilde.dml <- numeric(0)
  r.est <- numeric(0)
  for (k in 1:K) {
    # Define the indices for the k-th fold
    idx_start <- (k - 1) * fold_size + 1
    idx_end <- k * fold_size
    idx <- idx_start:idx_end
    
    # Select the test fold
    d.sel <- D[idx]
    x.sel <- X[idx,]
    y.sel <- Y[idx]
    
    if (est == "ols"){
      # Estimate on training data
      D_X <- as.vector(lm(D[-idx, ] ~ X[-idx, ])$coef)
      
      # Create prediction data (add intercept column)
      new <- cbind(1, x.sel)
      
      # Predict rhat
      rhat <- as.vector(D_X %*% t(new))  
    } else if (est == "nw"){
      # Split the data
      # Split the data
      fdat <- cbind(as.data.frame(D),as.data.frame(X))
      names(fdat) <- c("D", "X1", "X2", "X3")
      train <- fdat[-idx, ]   # new values to predict
      test <- fdat[idx, ]   # new values to predict
      
      # Fit N-W regression on training data
      nw_model <- npreg(D ~ X1 + X2 + X3, data = train, regtype = "lc", bwmethod = "cv.aic")
      
      # Predict on x.sel
      rhat <- predict(nw_model, newdata = test)
    } else if ( est == "poly" ){
      X <- as.matrix(X)
      d <- ncol(X)
      BX <- matrix(1, n, 1) 
      for (j in 1:d) {
        for (deg in 1:p) {
          BX <- cbind(BX, X[, j]^deg)
        }
      }
      D_X <- as.vector(lm(D[-idx] ~ -1+ BX[-idx, ])$coef)
      rhat <- BX[idx, ] %*% as.matrix(D_X) 
      
    }
    
    # saving estimated alpha values
    r.est <- c(r.est, rhat)
    
    # Construct dataset for dml estimation
    dml.dat <- data.frame(X = x.sel, D = d.sel, Y = y.sel, r = rhat, gx = gx[idx])
    
    # Apply the oracle DR estimator
    Psi_tilde.l.dml <- dml.est(dml.dat)
    
    # Store the result
    Psi_tilde.dml <- c(Psi_tilde.dml, Psi_tilde.l.dml)
  }
  
  theta.dml <- mean(Psi_tilde.dml)
  devs <- Psi_tilde.dml-theta.dml
  se.dml <- sqrt(mean((devs)^2)/n)
  
  return(list("est" =theta.dml, 
              "se" = se.dml, 
              "rest" = r.est))
}


dml.est <- function(dat) {
  gamma <- dat$gx        # true γ0(x)
  alpha <- 1 / dat$r           # true α0(x)
  score <- alpha * dat$D * (dat$Y - gamma) + gamma 
}

oracle_dr <- function(dat) {
  n <- nrow(dat$X)
  gamma <- dat$gx        # true γ0(x)
  alpha <- 1 / dat$r           # true α0(x)
  score <- alpha * dat$D * (dat$Y - gamma) + gamma 
  theta_oracle <- mean(score)
  devs <- score-theta_oracle
  se.oracle <- sqrt(mean((devs)^2)/n)
  return(c("est" =theta_oracle, 
           "se" = se.oracle))
}

generate_data <- function(n = 200, 
                          beta0 = 0.25, 
                          beta1 = 0.5, 
                          sigma = 1, 
                          model = "logit",
                          seed = 1) {
  set.seed(seed)
  
  # Define the true function g(x)
  g_fun <- function(x) exp(x)
  
  # Generate data
  X <- runif(n)               # Uniform[0, 1]
  V <- runif(n)
  U <- rnorm(n)              # N(0, 1)
  
  index <- beta0 + beta1 * X
  r <- if (model == "probit") pnorm(index) else 1 / (1 + exp(-index))
  D <- as.integer(r > V)
  gx <- g_fun(X)
  Y <- gx + sigma * U
  
  # Return data as a data.frame
  data.frame(X = X, D = D, Y = Y, r = r, gx = gx)
}


# Load required package for Toeplitz matrix
# For toeplitz function

generate_data_toeplitz <- function(n = 200, 
                                   p = 3, 
                                   a0 = 0.5, 
                                   a1 = 0.3, 
                                   s1 = 1, 
                                   alpha = 1, 
                                   b0 = 0.2, 
                                   b1 = 0.4, 
                                   s2 = 1, 
                                   model = "probit", 
                                   seed = 2162016) {
  set.seed(seed)
  
  # Create Toeplitz correlation matrix
  rho <- 0.7
  corr_matrix <- toeplitz(rho^(0:(p-1)))
  
  # Cholesky decomposition (upper triangular, like MATLAB's chol)
  S <- chol(corr_matrix)
  
  # Generate covariates x with correlation
  X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% S
  
  # Extract x1 and x3 for compatibility with previous DGP
  x1 <- X[, 1]
  x3 <- X[, 3]
  
  # Generate latent index for d
  index <- a0 * x1 + a1 * exp(x3) / (1 + exp(x3)) + s1 * rnorm(n)
  
  # Generate binary d using logit or probit
  r <- if (model == "probit") pnorm(index) else 1 / (1 + exp(-index))
  D <- as.integer(r > runif(n))
  
  # Generate y
  # y <- alpha * d 
  gx <- b0 * exp(x1) / (1 + exp(x1)) + b1 * x3
  Y <- gx + s2 * rnorm(n)
  
  # Return data as a data frame
  list(X = X, D = D, Y = Y, r = r, gx = gx)
}









