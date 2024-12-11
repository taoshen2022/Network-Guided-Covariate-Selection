source("./Functions/Functions-S2.R")
require(randomForest)

n = 1000; 
n1 = 1000;
n2 = 200;
repetition = 50;

L = 10; #the underlying dimension

# set up the covariate matrix
p = 1200;
signal_size = 50;
signals = sample(p, signal_size);
Xmean = matrix(0, nrow = n, ncol = p);

Time_record = matrix(0, nrow = repetition, ncol = 5) #store the true positives
colnames(Time_record) = c("NGfs", "NGfs2", "Marginal", "RF1", "Oracle")
for (ii in 1: repetition){
  mu = 0.3;
  
  # feature selection error
  corr <- matrix(0, nrow= L, ncol = L);
  corr[1:3,1:3] <- matrix(runif(9, 0,1.0), nrow = 3, ncol = 3);
  corr[4:6,4:6] <- matrix(runif(9, 0,1.0), nrow = 3, ncol = 3);
  corr[7:10,7:10] <- matrix(runif(16, 0,1), nrow = 4, ncol = 4);
  diag(corr) <- 1;
  corr <- as.matrix(forceSymmetric(corr));
  
  muvec <- c(rep(2,10));
  Y <- rmvnorm(n, mean = muvec, sigma = corr);
  
  Omega = Y%*%t(Y)*0.01;
  # Generate the covariates
  beta = matrix(abs(runif(L*signal_size, 0.05, mu)), 
                nrow = L, ncol = signal_size);
  beta =  beta * matrix(rbinom(L*signal_size, 1, 1/2)*2-1, nrow = L, ncol = signal_size);
  Xmean[,signals] = Y%*%beta;
  
  # set up the response variable
  gamma = matrix(abs(rnorm(L, 0, 1)), nrow = L, ncol = 1);
  gamma = gamma - mean(gamma);
  Zmean = Y%*%gamma;
  
  X = Xmean + matrix(rnorm(n*p), n, p); 
  Z = Zmean + rnorm(n, 0, sqrt(0.5));
  
  # Generate the network
  A = matrix(runif(n*n, 0, 1), nrow = n);
  A = Omega - A;
  A = 1*(A >= 0);
  diag(A) = 0;
  A <- as.matrix(forceSymmetric(A));
  A = A[1:n1, 1:n1]; #only the first n1 nodes are observed with the network
  
  
  ######### Covariate selection ###################
  
  # New approach: Network-Guided
  start_time <- Sys.time()
  fsresult = NGCS(A, X[1:n1, ], L, cri = cri, Lap = FALSE, Standardize = FALSE); 
  pval.new = fsresult$pval;
  end_time <- Sys.time()
  Time_record[ii, "NGfs"] = end_time - start_time;
  
  # New approach: Network-Guided, Laplacian
  start_time <- Sys.time()
  fslap = NGCS(A, X[1:n1, ], L, cri = cri, Lap = TRUE, Standardize = FALSE); 
  end_time <- Sys.time()
  pval.lap = fslap$pval;
  Time_record[ii, "NGfs2"] = end_time - start_time;
  
  # Marginal
  start_time <- Sys.time()
  pval.mar = apply(X, 2, function(x) regp(lm(Z[(n1+1):n]~x[(n1+1):n]-1)));
  end_time <- Sys.time()
  Time_record[ii, "Marginal"] = end_time - start_time;
  
  # random forest
  start_time <- Sys.time()
  data <- as.data.frame(cbind(Z,X));
  data.rf <- randomForest(V1 ~ ., data=data, ntree=100,
                          keep.forest=FALSE, importance=TRUE);
  rf.re <- t(data.rf$importance[,1]); 
  end_time <- Sys.time()
  Time_record[ii, "RF1"] = end_time - start_time;
  
  # Latent space known
  start_time <- Sys.time()
  pval.ls =  apply(X, 2, function(x) regp(lm(x[1:n1]~Y[1:n1,]-1)));
  end_time <- Sys.time()
  Time_record[ii, "Oracle"] = end_time - start_time;
  
  filename = paste('./Exp3_rep_', ii, '.Rdata', sep = "");
  save(pval.new, pval.lap, pval.mar, rf.re, pval.ls, file = filename)
  print(paste("rep = ", ii));
}

pseq <- seq(10,100,10);
TP2 = matrix(0, nrow = length(pseq), ncol = 5); #store the true positives
colnames(TP2) = c("NGfs", "NGfs2", "Marginal", "RF1", "Oracle");
P2 = TP2; #store the positives

for(mujj in 1:length(pseq)){
  cri = pseq[mujj]

  TPmat = matrix(0, nrow = repetition, ncol = 5); #store the #TP in each repetition
  colnames(TPmat) = c("NGfs", "NGfs2", "Marginal", "RF1", "Oracle");
  
  #Calculation of the error 
  for (ii in 1: repetition){
    #set.seed(ii+9999)
    
    sig_vec = rep(0, p); sig_vec[signals] = 1;
    
    load(paste('./Exp3_rep_', ii, '.Rdata', sep = ""))
    
    # New approach: Network-Guided
    fsresult = HCfs_cri(pval.new, cri = cri);
    TPmat[ii, "NGfs"] = sum(sig_vec[fsresult$selected]);
    
    # New approach: Network-Guided, Laplacian
    fslap = HCfs_cri(pval.lap, cri = cri);
    TPmat[ii, "NGfs2"] = sum(sig_vec[fslap$selected]);
    
    # Marginal
    marresult = HCfs_cri(pval.mar, cri = cri);
    TPmat[ii, "Marginal"] = sum(sig_vec[marresult$selected]);
    
    ##random forest
    rfresult1 <- which(rank(rf.re)>=p-cri+1);
    TPmat[ii, "RF1"] = sum(sig_vec[rfresult1]);
    
    # Latent space known
    oracle = HCfs_cri(pval.ls, cri = cri);
    TPmat[ii, "Oracle"] = sum(sig_vec[oracle$selected]);
    
    
    
    #print(ii)
  }
  
  TP2[mujj,] = colMeans(TPmat); 
  
  print(paste("p = ", cri));
}

pseq44 <-pseq
save(pseq44, TP2,file = "fig1-right-rdpg.Rdata")

