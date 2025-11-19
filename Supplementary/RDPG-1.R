source("./Functions/Functions-S2.R")

n = 1000; 
n1 = 800;
n2 = 200;
repetition = 50;

L = 10; #the underlying dimension

# set up the covariate matrix
p = 1200;
signal_size = 50;
signals = sample(p, signal_size);
Xmean = matrix(0, nrow = n, ncol = p);
mseq = seq(0.1, 1.0, 0.05);



tau_matrix = matrix(0, nrow = length(mseq), ncol = 5); #store the true positives
colnames(tau_matrix) = c("NGfs", "NGfs2", "NGfs2l", "NGfs22l", "Oracle");
mj_matrix <- tau_matrix;


#set up the difference between the different communities
for(ljj in 1:length(mseq)){
  
  L = 10;

  tau_mat = matrix(0, nrow = repetition, ncol = 5); #store the true positives
  colnames(tau_mat) = c("NGfs", "NGfs2", "NGfs2k", "NGfs22k", "Oracle");
  mj_mat <- tau_mat;
  time_mat <- tau_mat;
  
  for (ii in 1:repetition){
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
    M = matrix(0, nrow = L, ncol = p);
    mu = mseq[ljj];
    beta = matrix(runif(L*signal_size, 0.05, mu), 
                  nrow = L, ncol = signal_size);
    beta =  beta * matrix(rbinom(L*signal_size, 1, 1/2)*2-1, nrow = L, ncol = signal_size);
    M[,signals] = beta;
    Xmean[,signals] = Y%*%beta;
    X = Xmean + matrix(rnorm(n*p), n, p); 
    
    # set up the response variable
    gamma = matrix(abs(rnorm(L, 0, 1)), nrow = L, ncol = 1);
    Zmean = Y%*%gamma;
    Z = Zmean + rnorm(n, 0, sqrt(0.5));
    
    # Generate the network
    A = matrix(runif(n*n, 0, 1), nrow = n);
    A = Omega - A;
    A = 1*(A >= 0);
    diag(A) = 0;
    A <- as.matrix(forceSymmetric(A));
    
    start_time <- Sys.time()
    ##oracle case
    op_vec <- svd(Y, nu = L, nv = L);
    op_vec <- op_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(op_vec)%*%Y%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "Oracle"] <- tau;
    mj_mat[ii, "Oracle"] <- Mj;
    time_mat[ii, "Oracle"] <- end_time - start_time;
    
    start_time <- Sys.time()
    ad_vec <- svd(A, nu = L, nv = L);
    ad_vec <- ad_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(ad_vec)%*%Y%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "NGfs"] <- tau;
    mj_mat[ii, "NGfs"] <- Mj;
    time_mat[ii, "NGfs"] <- end_time - start_time;
    
    start_time <- Sys.time()
    ad_vec <- svd(A, nu = 2*L, nv = 2*L);
    ad_vec <- ad_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(ad_vec)%*%Y%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "NGfs2k"] <- tau;
    mj_mat[ii, "NGfs2k"] <- Mj;
    time_mat[ii, "NGfs2k"] <- end_time - start_time;
    
    #Laplacian case
    start_time <- Sys.time()
    Net = A;
    d = colSums(A); #degree vector
    ind_neighbor = which(d > 0); #indices for nodes with at least one 
    Dtau=d[ind_neighbor]^(-0.5);
    Ltau = t(Dtau * A[ind_neighbor, ind_neighbor]) * Dtau;
    Net[ind_neighbor, ind_neighbor] = Ltau;
    lp_vec <- svd(Net, nu = L, nv = L);
    lp_vec <- lp_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(lp_vec)%*%Y%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "NGfs2"] <- tau;
    mj_mat[ii, "NGfs2"] <- Mj;
    time_mat[ii, "NGfs2"] <- end_time - start_time;
    
    start_time <- Sys.time()
    lp_vec <- svd(Net, nu = 2*L, nv = 2*L);
    lp_vec <- lp_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(lp_vec)%*%Y%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "NGfs22k"] <- tau;
    mj_mat[ii, "NGfs22k"] <- Mj;
    time_mat[ii, "NGfs22k"] <- end_time - start_time;
    
  }
  filename = paste('./Exp1_mu_', mu, '.Rdata', sep = "");
  save(tau_mat, mj_mat, time_mat, file = filename)
  tau_matrix[ljj,] = colMeans(tau_mat); mj_matrix[ljj,] = colMeans(mj_mat)
  
  print(paste("L = ", L));
}

tau_matrix33 <- tau_matrix
save(mseq, tau_matrix33, file = "fig1-left-rgpg.Rdata")
