## Revision: verifying for effectiveness of the testing step


repetition <- 50
K          <- 3                    
a_seq      <- seq(0.1, 1, by = 0.2)

method_names <- c("NGfs-T", "NGfs2-T", "NGfs-H", "NGfs2-H",
                  "NGfs-G", "NGfs2-G", "Chi-square", "SFS")

FDR   <- matrix(NA_real_, nrow = length(a_seq), ncol = length(method_names))
HCtab <- matrix(0,        nrow = length(a_seq), ncol = length(method_names))
colnames(FDR)   <- method_names
colnames(HCtab) <- method_names

# HC threshold used for global detection
hc_thresh <- sqrt(2 * log(log(p)))

for (a_idx in seq_along(a_seq)) {
  a  <- a_seq[a_idx]
  mu <- 0.07          # fixed signal strength in this experiment
  
  FDRmat <- matrix(NA_real_, nrow = repetition, ncol = length(method_names))
  HCmat  <- matrix(0,        nrow = repetition, ncol = length(method_names))
  colnames(FDRmat) <- method_names
  colnames(HCmat)  <- method_names
  
  for (ii in seq_len(repetition)) {
    set.seed(1000 + ii)   
    
   
    de <- 0.05
    M  <- matrix(0, nrow = K, ncol = p)
    for (j in signals) {
      M[, j] <- rnorm(K, mean = mu, sd = de)
    }
    M <- M * matrix(rbinom(K * p, 1, 0.5) * 2 - 1, nrow = K, ncol = p)
    
    B <- matrix(c(1, a, a,
                  a, 1, a,
                  a, a, 1),
                ncol = 3, byrow = TRUE)
    
    l  <- sample(seq_len(K), n, replace = TRUE)  # node labels
    Pi <- matrix(0, nrow = n, ncol = K)
    for (k in seq_len(K)) {
      Pi[l == k, k] <- 1
    }
    
    theta <- rexp(n, rate = 5) + 0.06
    Theta <- diag(theta)
    
    Omega <- Theta %*% Pi %*% B %*% t(Pi) %*% Theta
    
    Xmean <- Pi %*% M
    X     <- Xmean + matrix(rnorm(n * p), nrow = n, ncol = p)
    
    A <- matrix(runif(n * n), nrow = n)
    A <- 1 * (Omega - A >= 0)
    diag(A) <- 0
    A <- as.matrix(forceSymmetric(A))
    A <- A[1:n1, 1:n1]    # only first n1 nodes observed with network
    
    cri     <- signal_size
    sig_vec <- integer(p); sig_vec[signals] <- 1
    
    # NGCS (adjacency)
    fs_adj <- NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE)
    pval.fs <- fs_adj$pval
    
    # NGCS (Laplacian)
    fs_lap <- NGCS(A, X[1:n1, ], K, Lap = TRUE, Standardize = FALSE)
    pval.fsl <- fs_lap$pval
    
    
    # (a) HCT on NGCS p-values
    fs_h <- HCfs(pval.fs)
    FDRmat[ii, "NGfs-H"] <-
      1 - sum(sig_vec[fs_h$selected]) / length(fs_h$selected)
    HCmat[ii, "NGfs-H"]  <-
      as.integer(max(fs_adj$hc) <= hc_thresh)
    
    # (b) HCT on NGCS(Lap) p-values
    fs_lap_h <- HCfs(pval.fsl)
    FDRmat[ii, "NGfs2-H"] <-
      1 - sum(sig_vec[fs_lap_h$selected]) / length(fs_lap_h$selected)
    HCmat[ii, "NGfs2-H"]  <-
      as.integer(max(fs_lap$hc) <= hc_thresh)
    
    # (c) Additional HCT
    fs_hp  <- HCfs(fs_adj$hpval)
    FDRmat[ii, "NGfs-G"] <-
      1 - sum(sig_vec[fs_hp$selected]) / length(fs_hp$selected)
    HCmat[ii, "NGfs-G"]  <-
      as.integer(max(fs_hp$HC) <= hc_thresh)
    
    
    cat("a =", a, ", rep =", ii, "\n")
  }
  

  for (m in seq_along(method_names)) {
    idx_detect <- which(HCmat[, m] == 0)   # runs where test rejects
    if (length(idx_detect) > 0) {
      FDR[a_idx, m] <- mean(FDRmat[idx_detect, m], na.rm = TRUE)
    } else {
      FDR[a_idx, m] <- NA_real_
    }
  }
  
  HCtab[a_idx, ] <- colSums(HCmat)
  
  cat("Finished a =", a, "\n")
}
