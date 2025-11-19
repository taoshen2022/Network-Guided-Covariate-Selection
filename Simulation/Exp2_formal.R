#Revised code: DCMM-based experiments

n  <- 1000
n1 <- 800      # number of observed nodes in the network

K  <- 3         # number of communities
p  <- 1200      # number of covariates

signal_size <- 50
signals     <- sample(p, signal_size)

museq      <- seq(0.1, 0.5, by = 0.05)
repetition <- 50

methods <- c("NGfs-T", "NGfs2-T",
             "NGfs-H", "NGfs2-H",
             "NGfs-G", "NGfs2-G",
             "Chi-square", "SFS")


## Local functions for assistance

make_B <- function(a = pi / 4) {
  matrix(
    c(
      sin(a)^2,      cos(a)^2 / 2, cos(a)^2 / 2,
      cos(a)^2 / 2,  sin(a)^2,     cos(a)^2 / 2,
      cos(a)^2 / 2,  cos(a)^2 / 2, sin(a)^2
    ),
    ncol = 3, byrow = TRUE
  )
}

generate_DCMM <- function(n, K, a = pi / 4,
                          rate = 5, shift = 0.06) {
  B <- make_B(a)
  
  # mixed-membership matrix Pi (DCMM)
  Pi <- matrix(0, n, K)
  for (i in seq_len(n)) {
    base_k <- sample(seq_len(K), 1)
    Pi[i, base_k] <- 1
  }
  Pi <- Pi + matrix(runif(n * K, 0, 0.3), n, K)
  Pi <- t(apply(Pi, 1, function(x) x / sum(x)))
  
  theta <- rexp(n, rate = rate) + shift
  Theta <- diag(theta)
  
  Omega <- Theta %*% Pi %*% B %*% t(Pi) %*% Theta
  
  list(Pi = Pi, Theta = Theta, Omega = Omega, B = B)
}

generate_M <- function(K, p, signals, mu, de = 0.05) {
  M <- matrix(0, nrow = K, ncol = p)
  for (j in signals) {
    M[, j] <- rnorm(K, mean = mu, sd = de)
  }
  # random signs
  M * matrix(rbinom(K * p, 1, 0.5) * 2 - 1, nrow = K, ncol = p)
}


## Case 1: DCMM + Gaussian noise


FDR_case1 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case1) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    dcmm  <- generate_DCMM(n, K)
    Pi    <- dcmm$Pi
    Omega <- dcmm$Omega
    
    M     <- generate_M(K, p, signals, mu)
    Xmean <- Pi %*% M
    
    X <- Xmean + matrix(rnorm(n * p), nrow = n, ncol = p)

    
    A <- matrix(runif(n * n), nrow = n)
    A <- 1 * (Omega - A >= 0)
    diag(A) <- 0
    A <- as.matrix(forceSymmetric(A))
    A <- A[1:n1, 1:n1]
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    # Network-guided screening (adjacency only)
    fs_adj <- NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE)
    pval.fs  <- fs_adj$pval
    
    # Network-guided screening with Laplacian
    fs_lap  <- NGCS(A, X[1:n1, ], K, Lap = TRUE, Standardize = FALSE)
    pval.fsl <- fs_lap$pval
    
    ## Fixed threshold
    fs_fix  <- HCfs_cri(pval.fs,  cri = cri)
    lap_fix <- HCfs_cri(pval.fsl, cri = cri)
    
    FDRmat[ii, "NGfs-T"]  <- 1 - sum(sig_vec[fs_fix$selected])  / length(fs_fix$selected)
    FDRmat[ii, "NGfs2-T"] <- 1 - sum(sig_vec[lap_fix$selected]) / length(lap_fix$selected)
    
    ## HCT (data-driven)
    fs_hc  <- HCfs(pval.fs)
    lap_hc <- HCfs(pval.fsl)
    
    FDRmat[ii, "NGfs-H"]  <- 1 - sum(sig_vec[fs_hc$selected])  / length(fs_hc$selected)
    FDRmat[ii, "NGfs2-H"] <- 1 - sum(sig_vec[lap_hc$selected]) / length(lap_hc$selected)
    
    ## HW-based p-values
    fs_hw  <- HCfs(fs_adj$hpval)
    lap_hw <- HCfs(fs_lap$hpval)
    
    FDRmat[ii, "NGfs-G"]  <- 1 - sum(sig_vec[fs_hw$selected])  / length(fs_hw$selected)
    FDRmat[ii, "NGfs2-G"] <- 1 - sum(sig_vec[lap_hw$selected]) / length(lap_hw$selected)
    
    ## Chi-square marginal
    pval.chi <- apply(X, 2, function(x) pchisq(sum(x^2), df = n, lower.tail = FALSE))
    chi_res  <- HCfs_cri(pval.chi, cri = cri)
    FDRmat[ii, "Chi-square"] <- 1 - sum(sig_vec[chi_res$selected]) / length(chi_res$selected)
    
    ## Sparse K-means feature selection
    skmeans   <- KMeansSparseCluster(X, K)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 1 (DCMM + Gaussian), mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case1[m_idx, ] <- colMeans(FDRmat)
  cat("Case 1 finished: mu =", mu, "\n")
}


## Case 2: DCMM + chi-square transform


FDR_case2 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case2) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    dcmm  <- generate_DCMM(n, K)
    Pi    <- dcmm$Pi
    Omega <- dcmm$Omega
    
    M     <- generate_M(K, p, signals, mu)
    Xmean <- Pi %*% M
    
    # chi-square sub-Gaussian transform 
    df <- 5
    E  <- matrix(rchisq(n * p, df = df), nrow = n, ncol = p)
    Z  <- ((E / df)^(1 / 3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))
    X  <- Xmean + Z
    
    A <- matrix(runif(n * n), nrow = n)
    A <- 1 * (Omega - A >= 0)
    diag(A) <- 0
    A <- as.matrix(forceSymmetric(A))
    A <- A[1:n1, 1:n1]
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    fs_adj <- NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE)
    pval.fs  <- fs_adj$pval
    
    fs_lap  <- NGCS(A, X[1:n1, ], K, Lap = TRUE, Standardize = FALSE)
    pval.fsl <- fs_lap$pval
    
    fs_fix  <- HCfs_cri(pval.fs,  cri = cri)
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
    
    skmeans   <- KMeansSparseCluster(X, K)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 2 (DCMM + chi-square), mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case2[m_idx, ] <- colMeans(FDRmat)
  cat("Case 2 finished: mu =", mu, "\n")
}


## Case 3: DCMM + sub-Gaussian (isotropic mixture)


FDR_case3 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case3) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    dcmm  <- generate_DCMM(n, K)
    Pi    <- dcmm$Pi
    Omega <- dcmm$Omega
    
    M     <- generate_M(K, p, signals, mu)
    Xmean <- Pi %*% M
    
    # Isotropic mixture noise
    noise <- rE_isotropic_mix(n = n, p = p, sigma = 1)$E
    X     <- Xmean + noise
    
    A <- matrix(runif(n * n), nrow = n)
    A <- 1 * (Omega - A >= 0)
    diag(A) <- 0
    A <- as.matrix(forceSymmetric(A))
    A <- A[1:n1, 1:n1]
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    fs_adj <- NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE)
    pval.fs  <- fs_adj$pval
    
    fs_lap  <- NGCS(A, X[1:n1, ], K, Lap = TRUE, Standardize = FALSE)
    pval.fsl <- fs_lap$pval
    
    fs_fix  <- HCfs_cri(pval.fs,  cri = cri)
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
    
    skmeans   <- KMeansSparseCluster(X, K)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 3 (DCMM + mix), mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case3[m_idx, ] <- colMeans(FDRmat)
  cat("Case 3 finished: mu =", mu, "\n")
}

