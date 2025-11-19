source("./Functions/Functions-S1.R")

repetition = 50;

K = 3; #number of communities

n = 1000; 
n1 = 800;
n2 = 200;
p = 1200;
signal_size = 50;
signals = sample(p, signal_size);

Time_record = matrix(0, nrow = repetition, ncol = 6) 
colnames(Time_record) = c("NGfs", "NGfs2", "Chi-square", "FSCA", "SFS", "Oracle")

for (ii in 1:repetition){
  mu = 0.25;
  de = 0.2;
  for (i in signals){
    M[,i] = matrix(runif(K, mu, mu+de), nrow = K, ncol = 1);
  }
  M = M*matrix(rbinom(K*p, 1, 1/2)*2-1, nrow = K, ncol = p);
  
  B = diag(rep(1-0.4, K)) + 0.4*ones(K);

  Pi = matrix(0, n, K); # label matrix
  for (i in 1:n){
    tar_index <- sample(1:3, 1);
    if (tar_index == 1){
      Pi[i, ] <- c(1,0,0);
    }
    if (tar_index == 2){
      Pi[i, ] <- c(0,1,0);
    }
    if (tar_index == 3){
      Pi[i, ] <- c(0,0,1);
    }
  }

  Pi <- Pi + matrix(runif(n*K, 0, 0.3), n, K);
  Pi <- t(apply(Pi, 1, function(x) x/sum(x)));
  
  theta = abs(rnorm(n, 0.2, sqrt(0.3)));
  Theta = diag(theta); # node degree heterogeneity
  
  # Get the expected adjacency matrix
  Omega = Theta %*% Pi %*% B %*% t(Pi) %*% Theta;
  
  # Expectation of the covariates
  Xmean = Pi%*%M;
  
  
  # Generate the covariates
  X = Xmean + matrix(rnorm(n*p), n, p); 
  
  # Generate the network
  A = matrix(runif(n*n, 0, 1), nrow = n);
  A = Omega - A;
  A = 1*(A >= 0);
  diag(A) = 0;
  A <- as.matrix(forceSymmetric(A));
  A = A[1:n1, 1:n1];
  
  
  ######### Covariate selection ###################
  
  # New approach: Network-Guided
  start_time <- Sys.time()
  fsresult = NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE);
  end_time <- Sys.time()
  pval.fs <- fsresult$pval;
  Time_record[ii,"NGfs"] <- end_time - start_time;
  
  # New approach: Network-Guided, Laplacian
  start_time <- Sys.time()
  fslap = NGCS(A, X[1:n1, ], K, Lap = TRUE, Standardize = FALSE); 
  end_time <- Sys.time()
  pval.fsl <- fslap$pval;
  Time_record[ii, "NGfs2"] <- end_time - start_time;
  
  # Chisquare
  start_time <- Sys.time()
  pval.chi = apply(X, 2, function(x) pchisq(sum(x^2), n, lower.tail = FALSE));
  end_time <- Sys.time()
  Time_record[ii, "Chi-square"] <- end_time - start_time;
  
  # FSCA
  start_time <- Sys.time()
  fsca.vx <- apply(X, 2, function(x) (1-norm(X-(x%*%t(x)/((t(x)%*%x)[1]))%*%X, type="F")^2/norm(X, type="F")^2)*100);
  end_time <- Sys.time()
  Time_record[ii, "FSCA"] <- end_time - start_time;
  
  # SFS
  start_time <- Sys.time()
  skmeans = KMeansSparseCluster(X, K);
  sk.vx = skmeans[[length(skmeans)]]$ws;
  end_time <- Sys.time()
  Time_record[ii, "SFS"] <- end_time - start_time;
  
  # Oracle
  start_time <- Sys.time()
  pval.label = apply(X[1:n1, ], 2, function(x)  regp(lm(x~as.factor(Pi[1:n1])-1)));
  end_time <- Sys.time()
  Time_record[ii, "Oracle"] <- end_time - start_time;
  
  
  filename = paste('./Exp3_rep_', ii, '.Rdata', sep = "");
  save(pval.fs, pval.fsl, pval.chi, fsca.vx, pval.label, sk.vx, file = filename)
  print(paste("rep = ", ii));
}


pseq = seq(10,100, by = 10); 

TP3 = matrix(0, nrow = length(pseq), ncol = 6) #store the true positives
colnames(TP3) = c("NGfs", "NGfs2", "Chi-square", "FSCA", "SFS", "Oracle")


#set up the difference between the different communities
for(pjj in 1:length(pseq)){
  cri = pseq[pjj];
  
  # feature selection error
  TPmat = matrix(0, nrow = repetition, ncol = 6); #store the #TP in each repetition
  colnames(TPmat) = c("NGfs", "NGfs2", "Chi-square", "FSCA", "SFS", "Oracle");
  
  
  for (ii in 1: repetition){
    load(paste('./Exp3_rep_', ii, '.Rdata', sep = ""));
    
    ######### Covariate selection ###################
    #Truth
    sig_vec = rep(0, p); sig_vec[signals] = 1;
    
    # New approach: Network-Guided
    fsresult = HCfs_cri(pval.fs, cri = cri);
    TPmat[ii, "NGfs"] = sum(sig_vec[fsresult$selected]);
    
    # New approach: Network-Guided, Laplacian
    fslap = HCfs_cri(pval.fsl, cri = cri);
    TPmat[ii, "NGfs2"] = sum(sig_vec[fslap$selected]);
    
    # Chisquare
    #pval.chi = apply(X, 2, function(x) pchisq(sum(x^2), n, lower.tail = FALSE));
    chiresult = HCfs_cri(pval.chi, cri = cri);
    TPmat[ii, "Chi-square"] = sum(sig_vec[chiresult$selected]);
    
    #FSCA
    fscaresult <- which(rank(fsca.vx)>=p-cri+1)
    TPmat[ii, "FSCA"] = sum(sig_vec[fscaresult]);
    
    # SFS
    skresult <- which(rank(sk.vx)>=p-cri+1)
    TPmat[ii, "SFS"] = sum(sig_vec[skresult]);
    
    # Oracle
    oracle = HCfs_cri(pval.label, cri = cri);
    TPmat[ii, "Oracle"] = sum(sig_vec[oracle$selected]);
    
    print(ii)
  }
  
  TP3[pjj,] = colMeans(TPmat);
  
  print(paste("p = ", cri));
}


TP33 <- TP3
pseq22 <- pseq
save(pseq22, TP33, file = "fig1-right-dcsbmm.Rdata")
