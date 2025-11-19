#Revised code: RDPG-based experiments

library(mvtnorm)

n  <- 1000
n1 <- 800           # number of nodes with observed network

L  <- 10             # latent RDPG dimension
p  <- 1200           # number of covariates

signal_size <- 50
signals     <- sample(p, signal_size)

museq      <- seq(0.05, 0.30, by = 0.05)
repetition <- 50

methods <- c("NGfs-T", "NGfs2-T",
             "NGfs-H", "NGfs2-H",
             "NGfs-G", "NGfs2-G",
             "Chi-square", "SFS")


## Local functions

make_corr <- function(L) {
  corr <- matrix(0, nrow = L, ncol = L)
  corr[1:3, 1:3]   <- matrix(runif(9,  0, 1.0), nrow = 3)
  corr[4:6, 4:6]   <- matrix(runif(9,  0, 1.0), nrow = 3)
  corr[7:10, 7:10] <- matrix(runif(16, 0, 1.0), nrow = 4)
  diag(corr) <- 1
  as.matrix(forceSymmetric(corr))
}


generate_RDPG <- function(n, L) {
  Sigma <- make_corr(L)
  muvec <- rep(0.2, L)
  Y     <- rmvnorm(n, mean = muvec, sigma = Sigma)
  Omega <- Y %*% t(Y) * 0.01
  list(Y = Y, Omega = Omega)
}


generate_XZ_mean <- function(Y, mu, signals, p, signal_size, L) {
  # coefficients for signal covariates
  beta <- matrix(
    abs(runif(L * signal_size, 0.02, mu)),
    nrow = L, ncol = signal_size
  )
  beta <- beta * matrix(rbinom(L * signal_size, 1, 0.5) * 2 - 1,
                        nrow = L, ncol = signal_size)
  
  Xmean <- matrix(0, nrow = nrow(Y), ncol = p)
  Xmean[, signals] <- Y %*% beta
  
  # response mean (not used in selection but part of the DGP)
  gamma <- matrix(abs(rnorm(L, 0, 1)), nrow = L, ncol = 1)
  gamma <- gamma - mean(gamma)
  Zmean <- drop(Y %*% gamma)
  
  list(Xmean = Xmean, Zmean = Zmean)
}


generate_A_from_Omega <- function(Omega, n1) {
  n <- nrow(Omega)
  A <- matrix(runif(n * n), nrow = n)
  A <- 1 * (Omega - A >= 0)
  diag(A) <- 0
  A <- as.matrix(forceSymmetric(A))
  A[1:n1, 1:n1]
}


## Case 1: RDPG + Gaussian noise


FDR_case1 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case1) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    rdpg <- generate_RDPG(n, L)
    Y     <- rdpg$Y
    Omega <- rdpg$Omega
    
    XZ    <- generate_XZ_mean(Y, mu, signals, p, signal_size, L)
    Xmean <- XZ$Xmean
    Zmean <- XZ$Zmean
    
    X <- Xmean + matrix(rnorm(n * p), nrow = n, ncol = p)
    Z <- Zmean + rnorm(n, mean = 0, sd = sqrt(0.5))  # not used further
    
    A <- generate_A_from_Omega(Omega, n1)

    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    # Network-guided screening (adjacency)
    fs_adj <- NGCS(A, X[1:n1, ], L, Lap = FALSE, Standardize = FALSE)
    pval.fs <- fs_adj$pval
    
    # Network-guided with Laplacian
    fs_lap  <- NGCS(A, X[1:n1, ], L, Lap = TRUE, Standardize = FALSE)
    pval.fsl <- fs_lap$pval
    
    # Fixed threshold 
    fs_fix  <- HCfs_cri(pval.fs, cri = cri)
    lap_fix <- HCfs_cri(pval.fsl, cri = cri)
    
    FDRmat[ii, "NGfs-T"]  <- 1 - sum(sig_vec[fs_fix$selected])  / length(fs_fix$selected)
    FDRmat[ii, "NGfs2-T"] <- 1 - sum(sig_vec[lap_fix$selected]) / length(lap_fix$selected)
    
    # HCT
    fs_hc  <- HCfs(pval.fs)
    lap_hc <- HCfs(pval.fsl)
    
    FDRmat[ii, "NGfs-H"]  <- 1 - sum(sig_vec[fs_hc$selected])  / length(fs_hc$selected)
    FDRmat[ii, "NGfs2-H"] <- 1 - sum(sig_vec[lap_hc$selected]) / length(lap_hc$selected)
    
    # HW-based p-values
    fs_hw  <- HCfs(fs_adj$hpval)
    lap_hw <- HCfs(fs_lap$hpval)
    
    FDRmat[ii, "NGfs-G"]  <- 1 - sum(sig_vec[fs_hw$selected])  / length(fs_hw$selected)
    FDRmat[ii, "NGfs2-G"] <- 1 - sum(sig_vec[lap_hw$selected]) / length(lap_hw$selected)
    
    # Chi-square marginal
    pval.chi <- apply(X, 2, function(x) pchisq(sum(x^2), df = n, lower.tail = FALSE))
    chi_res  <- HCfs_cri(pval.chi, cri = cri)
    FDRmat[ii, "Chi-square"] <- 1 - sum(sig_vec[chi_res$selected]) / length(chi_res$selected)
    
    # Sparse K-means feature selection
    skmeans   <- KMeansSparseCluster(X, L)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 1 (RDPG + Gaussian), mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case1[m_idx, ] <- colMeans(FDRmat)
  cat("Case 1 finished: mu =", mu, "\n")
}


## Case 2: RDPG + chi-square transform


FDR_case2 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case2) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    rdpg <- generate_RDPG(n, L)
    Y     <- rdpg$Y
    Omega <- rdpg$Omega
    
    XZ    <- generate_XZ_mean(Y, mu, signals, p, signal_size, L)
    Xmean <- XZ$Xmean
    Zmean <- XZ$Zmean
    
    # chi-square transform
    df <- 5
    E  <- matrix(rchisq(n * p, df = df), nrow = n, ncol = p)
    Znoise <- ((E / df)^(1 / 3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))
    
    X <- Xmean + Znoise
    Z <- Zmean + rnorm(n, mean = 0, sd = sqrt(0.5))
    
    A <- generate_A_from_Omega(Omega, n1)
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    fs_adj <- NGCS(A, X[1:n1, ], L, Lap = FALSE, Standardize = FALSE)
    pval.fs <- fs_adj$pval
    
    fs_lap  <- NGCS(A, X[1:n1, ], L, Lap = TRUE, Standardize = FALSE)
    pval.fsl <- fs_lap$pval
    
    fs_fix  <- HCfs_cri(pval.fs, cri = cri)
    lap_fix <- HCfs_cri(pval.fsl, cri = cri)
    
    FDRmat[ii, "NGfs-T"]  <- 1 - sum(sig_vec[fs_fix$selected])  / length(fs_fix$selected)
    FDRmat[ii, "NGfs2-T"] <- 1 - sum(sig_vec[lap_fix$selected]) / length(lap_fix$selected)
    
    fs_hc  <- HCfs(pval.fs)
    lap_hc <- HCfs(pval.fsl)
    
    FDRmat[ii, "NGfs-H"]  <- 1 - sum(sig_vec[fs_hc$selected])  / length(fs_hc$selected)
    FDRmat[ii, "NGfs2-H"] <- 1 - sum(sig_vec[lap_hc$selected]) / length(lap_hc$selected)
    
    fs_hw  <- HCfs(fs_adj$hpval)
    lap_hw <- HCfs(fs_lap$hpval)
    
    FDRmat[ii, "NGfs-G"]  <- 1 - sum(sig_vec[fs_hw$selected])  / length(fs_hw$selected)
    FDRmat[ii, "NGfs2-G"] <- 1 - sum(sig_vec[lap_hw$selected]) / length(lap_hw$selected)
    
    pval.chi <- apply(X, 2, function(x) pchisq(sum(x^2), df = n, lower.tail = FALSE))
    chi_res  <- HCfs_cri(pval.chi, cri = cri)
    FDRmat[ii, "Chi-square"] <- 1 - sum(sig_vec[chi_res$selected]) / length(chi_res$selected)
    
    skmeans   <- KMeansSparseCluster(X, L)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 2 (RDPG + chi-square), mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case2[m_idx, ] <- colMeans(FDRmat)
  cat("Case 2 finished: mu =", mu, "\n")
}



## Case 3: RDPG + sub-Gaussian (isotropic mixture)


FDR_case3 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case3) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    rdpg <- generate_RDPG(n, L)
    Y     <- rdpg$Y
    Omega <- rdpg$Omega
    
    XZ    <- generate_XZ_mean(Y, mu, signals, p, signal_size, L)
    Xmean <- XZ$Xmean
    Zmean <- XZ$Zmean
    
    # isotropic mixture noise
    noise <- rE_isotropic_mix(n = n, p = p, sigma = 1)$E
    X     <- Xmean + noise
    Z     <- Zmean + rnorm(n, mean = 0, sd = sqrt(0.5))
    
    A <- generate_A_from_Omega(Omega, n1)
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    fs_adj <- NGCS(A, X[1:n1, ], L, Lap = FALSE, Standardize = FALSE)
    pval.fs <- fs_adj$pval
    
    fs_lap  <- NGCS(A, X[1:n1, ], L, Lap = TRUE, Standardize = FALSE)
    pval.fsl <- fs_lap$pval
    
    fs_fix  <- HCfs_cri(pval.fs, cri = cri)
    lap_fix <- HCfs_cri(pval.fsl, cri = cri)
    
    FDRmat[ii, "NGfs-T"]  <- 1 - sum(sig_vec[fs_fix$selected])  / length(fs_fix$selected)
    FDRmat[ii, "NGfs2-T"] <- 1 - sum(sig_vec[lap_fix$selected]) / length(lap_fix$selected)
    
    fs_hc  <- HCfs(pval.fs)
    FDRmat[ii, "NGfs-H"] <- 1 - sum(sig_vec[fs_hc$selected]) / length(fs_hc$selected)
    
    lap_hc <- HCfs(pval.fsl)
    FDRmat[ii, "NGfs2-H"] <- 1 - sum(sig_vec[lap_hc$selected]) / length(lap_hc$selected)
    
    fs_hw  <- HCfs(fs_adj$hpval)
    lap_hw <- HCfs(fs_lap$hpval)
    
    FDRmat[ii, "NGfs-G"]  <- 1 - sum(sig_vec[fs_hw$selected])  / length(fs_hw$selected)
    FDRmat[ii, "NGfs2-G"] <- 1 - sum(sig_vec[lap_hw$selected]) / length(lap_hw$selected)
    
    # Chi-square marginal
    pval.chi <- apply(X, 2, function(x) pchisq(sum(x^2), df = n, lower.tail = FALSE))
    chi_res  <- HCfs_cri(pval.chi, cri = cri)
    FDRmat[ii, "Chi-square"] <- 1 - sum(sig_vec[chi_res$selected]) / length(chi_res$selected)
    
    # If you want SFS here as well, uncomment:
    # skmeans   <- KMeansSparseCluster(X, L)
    # sk_ws     <- skmeans[[length(skmeans)]]$ws
    # sk_select <- which(rank(sk_ws) >= p - cri + 1)
    # FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 3 (RDPG + mix), mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case3[m_idx, ] <- colMeans(FDRmat)
  cat("Case 3 finished: mu =", mu, "\n")
}


