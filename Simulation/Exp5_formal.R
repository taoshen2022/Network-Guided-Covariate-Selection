#Revised code: regression application (NG-reg implementation)


library(Matrix)
library(glmnet)
library(ncvreg)
library(care)     


n  <- 1000
n1 <- 800
n2 <- 200               # unused, kept for clarity
L  <- 10                # latent dimension for RDPG

p           <- 1200
signal_size <- 50
signals     <- sample(p, signal_size)    # fixed signal locations across all cases

mseq       <- c(0.5, 1.0, 2.0)          # signal strength for beta
repetition <- 50

method_names <- c("Lasso", "NGreg", "NGreg 2",
                  "SLM", "PCR", "MCP", "SCAD", "Base")

idx_test  <- 1:n1
idx_train <- (n1 + 1):n
n_train   <- length(idx_train)

## Local functions

make_corr <- function(L) {
  corr <- matrix(0, nrow = L, ncol = L)
  corr[1:3, 1:3]   <- matrix(runif(9,  0, 1.0), nrow = 3)
  corr[4:6, 4:6]   <- matrix(runif(9,  0, 1.0), nrow = 3)
  corr[7:10, 7:10] <- matrix(runif(16, 0, 1.0), nrow = 4)
  diag(corr) <- 1
  as.matrix(forceSymmetric(corr))
}


simulate_RDPG_regression <- function(mu_signal,
                                     noise_type = c("gaussian", "chisq", "isotropic"),
                                     mu_latent = 0.2) {
  noise_type <- match.arg(noise_type)
  
  ## latent positions 
  Sigma_Y <- make_corr(L)
  mu_vec  <- rep(mu_latent, L)
  Y       <- mvtnorm::rmvnorm(n, mean = mu_vec, sigma = Sigma_Y)
  
  Omega <- Y %*% t(Y) * 0.01
  A     <- matrix(runif(n * n), nrow = n)
  A     <- 1 * (Omega - A >= 0)
  diag(A) <- 0
  A <- as.matrix(forceSymmetric(A))
  A <- A[idx_test, idx_test]     
  
  Xmean <- matrix(0, nrow = n, ncol = p)
  beta  <- matrix(abs(runif(L * signal_size, 0.05, mu_signal)),
                  nrow = L, ncol = signal_size)
  beta  <- beta * matrix(rbinom(L * signal_size, 1, 0.5) * 2 - 1,
                         nrow = L, ncol = signal_size)
  Xmean[, signals] <- Y %*% beta
  
  ## covariate noise
  if (noise_type == "gaussian") {
    X <- Xmean + matrix(rnorm(n * p, 0, sqrt(0.5)), nrow = n, ncol = p)
  } else if (noise_type == "chisq") {
    df <- 5
    E  <- matrix(rchisq(n * p, df = df), nrow = n, ncol = p)
    Zt <- ((E / df)^(1 / 3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))
    X  <- Xmean + Zt
  } else if (noise_type == "isotropic") {
    X <- Xmean + r_isotropic(n = n, p = p, sigma = 1)
  }
  
  ## response
  gamma  <- matrix(abs(rnorm(L, 0, 1)), nrow = L)
  gamma  <- gamma - mean(gamma)
  Zmean  <- as.vector(Y %*% gamma)
  Z      <- Zmean + rnorm(n, 0, sqrt(0.5))
  
  list(A = A, X = X, Z = Z, Zmean = Zmean)
}



run_regression_once <- function(A, X, Z, Zmean) {
  out_mse  <- setNames(numeric(length(method_names)), method_names)
  out_mse2 <- out_mse
  
  fs1 <- NGCS(A, X[idx_test, ], L,   Lap = FALSE, Standardize = FALSE)
  fs2 <- NGCS(A, X[idx_test, ], 2*L, Lap = FALSE, Standardize = FALSE)
  
  ## --- Lasso  ---
  cv_lasso <- cv.glmnet(X[idx_train, ], Z[idx_train],
                        alpha = 1, family = "gaussian")
  est_lasso <- as.numeric(predict(cv_lasso,
                                  newx = X[idx_test, ], s = "lambda.min"))
  out_mse["Lasso"]  <- mean((est_lasso - Z[idx_test])^2)
  out_mse2["Lasso"] <- mean((est_lasso - Zmean[idx_test])^2)
  
  ## --- NGreg  ---
  est_fs <- reg_pred(X[idx_train, fs1$selected, drop = FALSE],
                     X[idx_test,  fs1$selected, drop = FALSE],
                     Z[idx_train],
                     k = min(L, fs1$S))
  out_mse["NGreg"]  <- mean((est_fs - Z[idx_test])^2)
  out_mse2["NGreg"] <- mean((est_fs - Zmean[idx_test])^2)
  
  ## --- NGreg (2L) ---
  est_fs2 <- reg_pred(X[idx_train, fs2$selected, drop = FALSE],
                      X[idx_test,  fs2$selected, drop = FALSE],
                      Z[idx_train],
                      k = min(L, fs2$S))
  out_mse["NGreg 2"]  <- mean((est_fs2 - Z[idx_test])^2)
  out_mse2["NGreg 2"] <- mean((est_fs2 - Zmean[idx_test])^2)
  
  ## --- PCR  ---
  est_pcr <- reg_pred(X[idx_train, ], X[idx_test, ], Z[idx_train], k = L)
  out_mse["PCR"]  <- mean((est_pcr - Z[idx_test])^2)
  out_mse2["PCR"] <- mean((est_pcr - Zmean[idx_test])^2)
  
  ## --- CAR + SLM ---
  car      <- carscore(X[idx_train, ], Z[idx_train], lambda = 0)
  car_pval <- 1 - pbeta(car^2, shape1 = 1/2,
                        shape2 = (n_train - 2) / 2)  # df based on training size
  car_sel  <- HCfs(car_pval)
  car_fit  <- slm(X[idx_train, car_sel$selected, drop = FALSE],
                  Z[idx_train], lambda = 0, lambda.var = 0)
  est_car  <- as.numeric(predict(car_fit,
                                 as.matrix(X[idx_test, car_sel$selected, drop = FALSE])))
  out_mse["SLM"]  <- mean((est_car - Z[idx_test])^2)
  out_mse2["SLM"] <- mean((est_car - Zmean[idx_test])^2, na.rm = TRUE)
  
  ## --- MCP / SCAD  ---
  mcp_model  <- cv.ncvreg(X[idx_train, ], Z[idx_train],
                          family = "gaussian", penalty = "MCP")
  scad_model <- cv.ncvreg(X[idx_train, ], Z[idx_train],
                          family = "gaussian", penalty = "SCAD")
  
  mcp_pred  <- as.numeric(predict(mcp_model,  X[idx_test, ], s = "lambda.min"))
  scad_pred <- as.numeric(predict(scad_model, X[idx_test, ], s = "lambda.min"))
  
  out_mse["MCP"]   <- mean((mcp_pred  - Z[idx_test])^2)
  out_mse2["MCP"]  <- mean((mcp_pred  - Zmean[idx_test])^2,  na.rm = TRUE)
  out_mse["SCAD"]  <- mean((scad_pred - Z[idx_test])^2)
  out_mse2["SCAD"] <- mean((scad_pred - Zmean[idx_test])^2, na.rm = TRUE)
  
  ## --- Using mean for prediction ---
  base_mean <- mean(Z[idx_train])
  base_pred <- rep(base_mean, length(idx_test))
  out_mse["Base"]  <- mean((base_pred - Z[idx_test])^2)
  out_mse2["Base"] <- mean((base_pred - Zmean[idx_test])^2)
  
  list(mse = out_mse, mse2 = out_mse2)
}



run_regression_case <- function(noise_type, mu_latent) {
  mse_mat  <- matrix(0, nrow = length(mseq), ncol = length(method_names))
  mse2_mat <- mse_mat
  colnames(mse_mat)  <- method_names
  colnames(mse2_mat) <- method_names
  
  for (m_idx in seq_along(mseq)) {
    mu_signal <- mseq[m_idx]
    
    msemat_rep  <- matrix(0, nrow = repetition, ncol = length(method_names))
    msemat2_rep <- msemat_rep
    colnames(msemat_rep)  <- method_names
    colnames(msemat2_rep) <- method_names
    
    for (ii in seq_len(repetition)) {
      set.seed(ii + 9999)
      
      sim <- simulate_RDPG_regression(mu_signal = mu_signal,
                                      noise_type = noise_type,
                                      mu_latent  = mu_latent)
      
      res <- run_regression_once(sim$A, sim$X, sim$Z, sim$Zmean)
      msemat_rep[ii, ]  <- res$mse
      msemat2_rep[ii, ] <- res$mse2
      
      cat("Noise =", noise_type,
          ", mu_signal =", mu_signal,
          ", rep =", ii, "\n")
    }
    
    mse_mat[m_idx, ]  <- colMeans(msemat_rep)
    mse2_mat[m_idx, ] <- colMeans(msemat2_rep)
    
    cat("Finished mu_signal =", mu_signal,
        "for noise =", noise_type, "\n")
  }
  
  list(mse = mse_mat, mse2 = mse2_mat)
}


## Run the three cases
res_reg_case1 <- run_regression_case("gaussian", mu_latent = 0.2)
res_reg_case2 <- run_regression_case("chisq", mu_latent = 0.2)
res_reg_case3 <- run_regression_case("isotropic", mu_latent = 0.2)

