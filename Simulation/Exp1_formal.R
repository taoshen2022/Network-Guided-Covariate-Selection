#Revised code: DCSBM-based experiments

n        <- 1000
n1       <- 800       # number of observed nodes in the network
K        <- 3          # number of communities
p        <- 1200       # number of covariates

signal_size <- 50
signals     <- sample(p, signal_size)

museq      <- seq(0.1, 0.5, by = 0.05)
repetition <- 50

methods <- c("NGfs-T", "NGfs2-T",
             "NGfs-H", "NGfs2-H",
             "NGfs-G", "NGfs2-G",
             "Chi-square", "SFS")

##Local functions for assistance
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

generate_DCSBM <- function(n, K, a = pi / 4, rate = 5, shift = 0.06) {
  B <- make_B(a)
  
  labels <- sample(seq_len(K), n, replace = TRUE)
  
  Pi <- matrix(0, n, K)
  for (k in seq_len(K)) {
    Pi[labels == k, k] <- 1
  }
  
  theta <- rexp(n, rate = rate) + shift
  Theta <- diag(theta)
  
  Omega <- Theta %*% Pi %*% B %*% t(Pi) %*% Theta
  
  A <- matrix(runif(n * n), nrow = n)
  A <- 1 * (Omega - A >= 0)
  diag(A) <- 0
  A <- as.matrix(forceSymmetric(A))
  
  list(A = A, Pi = Pi, Theta = Theta, Omega = Omega, labels = labels)
}

generate_M <- function(K, p, signals, mu, de = 0.05) {
  M <- matrix(0, nrow = K, ncol = p)
  for (j in signals) {
    M[, j] <- rnorm(K, mean = mu, sd = de)
  }

  M * matrix(rbinom(K * p, 1, 0.5) * 2 - 1, nrow = K, ncol = p)
}



## Case 1: DCSBM + Gaussian noise

FDR_case1 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case1) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    dcsbm <- generate_DCSBM(n, K)
    A     <- dcsbm$A[1:n1, 1:n1]     # only first n1 nodes observed
    Pi    <- dcsbm$Pi
    
    M     <- generate_M(K, p, signals, mu)
    Xmean <- Pi %*% M
    X     <- Xmean + matrix(rnorm(n * p), nrow = n, ncol = p)
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    ## NGCS (adjacency + Laplacian)
    fs_adj <- NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE)
    fs_lap <- NGCS(A, X[1:n1, ], K, Lap = TRUE,  Standardize = FALSE)
    
    pval.fs  <- fs_adj$pval
    pval.fsl <- fs_lap$pval
    
    ## Fixed threshold (top-cri via HCfs_cri)
    fs_fix   <- HCfs_cri(pval.fs,  cri = cri)
    lap_fix  <- HCfs_cri(pval.fsl, cri = cri)
    
    FDRmat[ii, "NGfs-T"]  <- 1 - sum(sig_vec[fs_fix$selected])  / length(fs_fix$selected)
    FDRmat[ii, "NGfs2-T"] <- 1 - sum(sig_vec[lap_fix$selected]) / length(lap_fix$selected)
    
    ## Higher Criticism thresholding 
    fs_hc   <- HCfs(pval.fs)
    lap_hc  <- HCfs(pval.fsl)
    
    FDRmat[ii, "NGfs-H"]  <- 1 - sum(sig_vec[fs_hc$selected])   / length(fs_hc$selected)
    FDRmat[ii, "NGfs2-H"] <- 1 - sum(sig_vec[lap_hc$selected])  / length(lap_hc$selected)
    
    ## HW-based p-values
    fs_hw   <- HCfs(fs_adj$hpval)
    lap_hw  <- HCfs(fs_lap$hpval)
    
    FDRmat[ii, "NGfs-G"]  <- 1 - sum(sig_vec[fs_hw$selected])   / length(fs_hw$selected)
    FDRmat[ii, "NGfs2-G"] <- 1 - sum(sig_vec[lap_hw$selected])  / length(lap_hw$selected)
    
    ## Chi-square marginal screening
    pval.chi  <- apply(X, 2, function(x) pchisq(sum(x^2), df = n, lower.tail = FALSE))
    chi_res   <- HCfs_cri(pval.chi, cri = cri)
    FDRmat[ii, "Chi-square"] <- 1 - sum(sig_vec[chi_res$selected]) / length(chi_res$selected)
    
    ## Ranking based on sparse K-means
    skmeans   <- KMeansSparseCluster(X, K)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 1, mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case1[m_idx, ] <- colMeans(FDRmat)
  cat("Case 1 finished: mu =", mu, "\n")
}


## Case 2: DCSBM + chi-square transform


FDR_case2 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case2) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    dcsbm <- generate_DCSBM(n, K)
    A     <- dcsbm$A[1:n1, 1:n1]
    Pi    <- dcsbm$Pi
    
    M     <- generate_M(K, p, signals, mu)
    Xmean <- Pi %*% M
    
    ## Sub-Gaussian transform of chi-square
    df    <- 10
    E     <- matrix(rchisq(n * p, df = df), nrow = n, ncol = p)
    Z     <- ((E / df)^(1 / 3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))
    X     <- Xmean + Z
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    fs_adj <- NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE)
    fs_lap <- NGCS(A, X[1:n1, ], K, Lap = TRUE,  Standardize = FALSE)
    
    pval.fs  <- fs_adj$pval
    pval.fsl <- fs_lap$pval
    
    fs_fix   <- HCfs_cri(pval.fs,  cri = cri)
    lap_fix  <- HCfs_cri(pval.fsl, cri = cri)
    
    FDRmat[ii, "NGfs-T"]  <- 1 - sum(sig_vec[fs_fix$selected])  / length(fs_fix$selected)
    FDRmat[ii, "NGfs2-T"] <- 1 - sum(sig_vec[lap_fix$selected]) / length(lap_fix$selected)
    
    fs_hc   <- HCfs(pval.fs)
    lap_hc  <- HCfs(pval.fsl)
    
    FDRmat[ii, "NGfs-H"]  <- 1 - sum(sig_vec[fs_hc$selected])   / length(fs_hc$selected)
    FDRmat[ii, "NGfs2-H"] <- 1 - sum(sig_vec[lap_hc$selected])  / length(lap_hc$selected)
    
    fs_hw   <- HCfs(fs_adj$hpval)
    lap_hw  <- HCfs(fs_lap$hpval)
    
    FDRmat[ii, "NGfs-G"]  <- 1 - sum(sig_vec[fs_hw$selected])   / length(fs_hw$selected)
    FDRmat[ii, "NGfs2-G"] <- 1 - sum(sig_vec[lap_hw$selected])  / length(lap_hw$selected)
    
    pval.chi <- apply(X, 2, function(x) pchisq(sum(x^2), df = n, lower.tail = FALSE))
    chi_res  <- HCfs_cri(pval.chi, cri = cri)
    FDRmat[ii, "Chi-square"] <- 1 - sum(sig_vec[chi_res$selected]) / length(chi_res$selected)
    
    skmeans   <- KMeansSparseCluster(X, K)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 2, mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case2[m_idx, ] <- colMeans(FDRmat)
  cat("Case 2 finished: mu =", mu, "\n")
}


## Case 3: DCSBM + sub-Gaussian (isotropic mixture)


FDR_case3 <- matrix(0, nrow = length(museq), ncol = length(methods))
colnames(FDR_case3) <- methods

for (m_idx in seq_along(museq)) {
  mu <- museq[m_idx]
  
  FDRmat <- matrix(0, nrow = repetition, ncol = length(methods))
  colnames(FDRmat) <- methods
  
  for (ii in seq_len(repetition)) {
    
    dcsbm <- generate_DCSBM(n, K)
    A     <- dcsbm$A[1:n1, 1:n1]
    Pi    <- dcsbm$Pi
    
    M     <- generate_M(K, p, signals, mu)
    Xmean <- Pi %*% M
    
    ## Sub-Gaussian isotropic mixture noise
    noise <- rE_isotropic_mix(n = n, p = p, sigma = 1)$E
    X     <- Xmean + noise
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    fs_adj <- NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE)
    fs_lap <- NGCS(A, X[1:n1, ], K, Lap = TRUE,  Standardize = FALSE)
    
    pval.fs  <- fs_adj$pval
    pval.fsl <- fs_lap$pval
    
    fs_fix   <- HCfs_cri(pval.fs,  cri = cri)
    lap_fix  <- HCfs_cri(pval.fsl, cri = cri)
    
    FDRmat[ii, "NGfs-T"]  <- 1 - sum(sig_vec[fs_fix$selected])  / length(fs_fix$selected)
    FDRmat[ii, "NGfs2-T"] <- 1 - sum(sig_vec[lap_fix$selected]) / length(lap_fix$selected)
    
    fs_hc   <- HCfs(pval.fs)
    lap_hc  <- HCfs(pval.fsl)
    
    FDRmat[ii, "NGfs-H"]  <- 1 - sum(sig_vec[fs_hc$selected])   / length(fs_hc$selected)
    FDRmat[ii, "NGfs2-H"] <- 1 - sum(sig_vec[lap_hc$selected])  / length(lap_hc$selected)
    
    fs_hw   <- HCfs(fs_adj$hpval)
    lap_hw  <- HCfs(fs_lap$hpval)
    
    FDRmat[ii, "NGfs-G"]  <- 1 - sum(sig_vec[fs_hw$selected])   / length(fs_hw$selected)
    FDRmat[ii, "NGfs2-G"] <- 1 - sum(sig_vec[lap_hw$selected])  / length(lap_hw$selected)
    
    pval.chi <- apply(X, 2, function(x) pchisq(sum(x^2), df = n, lower.tail = FALSE))
    chi_res  <- HCfs_cri(pval.chi, cri = cri)
    FDRmat[ii, "Chi-square"] <- 1 - sum(sig_vec[chi_res$selected]) / length(chi_res$selected)
    
    skmeans   <- KMeansSparseCluster(X, K)
    sk_ws     <- skmeans[[length(skmeans)]]$ws
    sk_select <- which(rank(sk_ws) >= p - cri + 1)
    
    FDRmat[ii, "SFS"] <- 1 - sum(sig_vec[sk_select]) / length(sk_select)
    
    cat("Case 3, mu =", mu, ", rep =", ii, "\n")
  }
  
  FDR_case3[m_idx, ] <- colMeans(FDRmat)
  cat("Case 3 finished: mu =", mu, "\n")
}


