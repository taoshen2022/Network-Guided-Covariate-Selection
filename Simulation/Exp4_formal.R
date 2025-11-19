#Revised code: clustering application (NG-clu implementation)


library(Matrix)
library(KMeansSparseCluster)

n  <- 1000
n1 <- 800
n2 <- 200            # currently unused, kept for completeness

K  <- 3              # number of communities
p  <- 1200
signal_size <- 50

museq      <- c(0.1, 0.3, 0.5)
repetition <- 50

method_names <- c("All covariates",
                  "NG clustering",
                  "NG clustering 2",
                  "IF-PCA",
                  "SKmeans",
                  "SAS")

## Local functions

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


simulate_DCSBM_cov <- function(n, n1, K, p,
                               signals,
                               mu, de,
                               noise_type = c("gaussian", "chisq", "isotropic")) {
  noise_type <- match.arg(noise_type)
  
  M <- matrix(0, nrow = K, ncol = p)
  for (j in signals) {
    M[, j] <- rnorm(K, mean = mu, sd = de)
  }
  M <- M * matrix(rbinom(K * p, 1, 0.5) * 2 - 1, nrow = K, ncol = p)
  
  B <- make_B()
  
  # labels and membership
  l <- sample(seq_len(K), n, replace = TRUE)
  Pi <- matrix(0, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    Pi[l == k, k] <- 1
  }
  
  # degree heterogeneity
  theta <- abs(rnorm(n, mean = 0.1, sd = sqrt(0.2)))
  Theta <- diag(theta)
  
  Omega <- Theta %*% Pi %*% B %*% t(Pi) %*% Theta
  
  Xmean <- Pi %*% M
  
  # covariate noise
  if (noise_type == "gaussian") {
    X <- Xmean + matrix(rnorm(n * p), nrow = n, ncol = p)
  } else if (noise_type == "chisq") {
    df <- 5
    E  <- matrix(rchisq(n * p, df = df), nrow = n, ncol = p)
    Z  <- ((E / df)^(1 / 3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))
    X  <- Xmean + Z
  } else if (noise_type == "isotropic") {
    
    X <- Xmean + r_isotropic(n = n, p = p, sigma = 1)
  }
  
  A <- matrix(runif(n * n), nrow = n)
  A <- 1 * (Omega - A >= 0)
  diag(A) <- 0
  A <- as.matrix(forceSymmetric(A))
  A <- A[1:n1, 1:n1]
  
  list(A = A, X = X, l = l)
}



run_clustering_once <- function(A, X, l, K) {
  p <- ncol(X)
  out_err <- setNames(numeric(length(method_names)), method_names)
  out_nmi <- out_err
  
  # NGCS (adjacency)
  fs_adj <- NGCS(A, X[1:nrow(A), ], K,   Lap = FALSE, Standardize = FALSE)
  # NGCS (2K para)
  fs_2K  <- NGCS(A, X[1:nrow(A), ], 2*K, Lap = FALSE, Standardize = FALSE)
  
  # chi-square
  pval.chi <- apply(X, 2, function(x) pchisq(sum(x^2), df = nrow(X), lower.tail = FALSE))
  chi_sel <- HCfs_cluster(pval.chi, nrow(X))
  
  ## ---- Clustering: All covariates ----
  est_all <- PCAclu(X, K)
  out_err["All covariates"] <- cluster(table(est_all, l))$error
  out_nmi["All covariates"] <- NMI(est_all, l)
  
  ## ---- Clustering: NGCS ----
  est_ng <- PCAclu(X[, fs_adj$selected, drop = FALSE], K)
  out_err["NG clustering"] <- cluster(table(est_ng, l))$error
  out_nmi["NG clustering"] <- NMI(est_ng, l)
  
  ## ---- Clustering: NGCS (2K) ----
  est_ng2 <- PCAclu(X[, fs_2K$selected, drop = FALSE], K)
  out_err["NG clustering 2"] <- cluster(table(est_ng2, l))$error
  out_nmi["NG clustering 2"] <- NMI(est_ng2, l)
  
  ## ---- Clustering: IF-PCA ----
  est_if <- PCAclu(X[, chi_sel$selected, drop = FALSE], K)
  out_err["IF-PCA"] <- cluster(table(est_if, l))$error
  out_nmi["IF-PCA"] <- NMI(est_if, l)
  
  ## ---- Clustering: sparse k-means ----
  sk <- KMeansSparseCluster(X, K)
  sk_cl <- sk[[length(sk)]]$Cs
  out_err["SKmeans"] <- cluster(table(sk_cl, l))$error
  out_nmi["SKmeans"] <- NMI(sk_cl, l)
  
  ## ---- Clustering: SAS ----
  center0 <- colMeans(X)
  Xc0     <- t(apply(X, 1, function(x) x - center0))
  tot     <- apply(Xc0, 2, function(x) sum(x^2))
  
  wcss <- numeric(p)
  for (j in seq_len(p)) {
    kmj <- kmeans(X[, j], centers = K)
    wcss[j] <- kmj$tot.withinss
  }
  score   <- (tot - wcss) / tot
  s       <- 50
  rank0   <- rank(score, ties.method = "random")
  init_set <- which(rank0 > p - s)
  
  sas_out <- sas(X, k = K, tot = tot, initial_set = init_set, s = s, itermax = 10)
  sas_cl  <- sas_out$result
  
  out_err["SAS"] <- cluster(table(sas_cl, l))$error
  out_nmi["SAS"] <- NMI(sas_cl, l)
  
  list(error = out_err, nmi = out_nmi)
}


run_clustering_case <- function(noise_type) {
  error_mat <- matrix(0, nrow = length(museq), ncol = length(method_names))
  colnames(error_mat) <- method_names
  nmi_mat <- error_mat
  
  for (m_idx in seq_along(museq)) {
    mu <- museq[m_idx]
    signals <- sample(p, signal_size)
    de <- 0.2
    
    err_rep <- matrix(0, nrow = repetition, ncol = length(method_names))
    colnames(err_rep) <- method_names
    nmi_rep <- err_rep
    
    for (ii in seq_len(repetition)) {
      set.seed(ii + 9999)
      
      sim <- simulate_DCSBM_cov(
        n      = n,
        n1     = n1,
        K      = K,
        p      = p,
        signals = signals,
        mu     = mu,
        de     = de,
        noise_type = noise_type
      )
      
      res <- run_clustering_once(sim$A, sim$X, sim$l, K)
      err_rep[ii, ] <- res$error
      nmi_rep[ii, ] <- res$nmi
      
      cat("Noise =", noise_type, ", mu =", mu, ", rep =", ii, "\n")
    }
    
    error_mat[m_idx, ] <- colMeans(err_rep)
    nmi_mat[m_idx, ]   <- apply(nmi_rep, 2, mean)
    
    cat("Finished mu =", mu, "for noise =", noise_type, "\n")
  }
  
  list(error = error_mat, nmi = nmi_mat)
}


## Run the three cases


res_case1 <- run_clustering_case("gaussian")
res_case2 <- run_clustering_case("chisq")
res_case3 <- run_clustering_case("isotropic")

# save(museq,
#      res_case1, res_case2, res_case3,
#      file = "clustering_DCSBM_results.RData")
