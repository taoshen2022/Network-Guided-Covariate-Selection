source("./Functions/Functions-S1.R")

repetition = 50;

K = 3; #number of communities
aseq = seq(0,2*pi,pi/16)

n = 1000; 
n1 = 800;
n2 = 200;
signal_size = 50;
signals = sample(p, signal_size);

tau_matrix = matrix(0, nrow = length(aseq), ncol = 5); #store the true positives
colnames(tau_matrix) = c("NGfs", "NGfs2", "NGfs2k", "NGfs22k", "Oracle");
mj_matrix <- tau_matrix;

#set up the difference between the different communities
for(ajj in 1:length(aseq)){
  
  
  a = aseq[ajj];
  mu = 0.30;
  de = 0.05;

  
  tau_mat = matrix(0, nrow = repetition, ncol = 5); #store the true positives
  colnames(tau_mat) = c("NGfs", "NGfs2", "NGfs2k", "NGfs22k", "Oracle");
  
  mj_mat <- tau_mat;
  time_mat <- tau_mat;
  
  for (ii in 1:repetition){
    M = matrix(0, nrow = K, ncol = p);
    for (i in signals){
      M[,i] = matrix(rnorm(K, mu, de), nrow = K, ncol = 1);
    }
    M = M*matrix(rbinom(K*p, 1, 1/2)*2-1, nrow = K, ncol = p);
    
    
    B = matrix(c(sin(a)^2, cos(a)^2/2, cos(a)^2/2, cos(a)^2/2, sin(a)^2, cos(a)^2/2, cos(a)^2/2, cos(a)^2/2, sin(a)^2), ncol=3, byrow = TRUE);
    
    
    l = sample(1:K, n, replace=TRUE); # node labels
    
    Pi = matrix(0, n, K); # label matrix
    for (k in 1:K){
      Pi[l == k, k] = 1;
    }
    
    theta = abs(rnorm(n, 0.1, sqrt(0.2)));

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
    
    start_time <- Sys.time()
    ##optimal case
    op_vec <- svd(Pi, nu = K, nv = K);
    op_vec <- op_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(op_vec)%*%Pi%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "Oracle"] <- tau;
    mj_mat[ii, "Oracle"] <- Mj;
    time_mat[ii, "Oracle"] <- end_time - start_time;
    
    
    ##adjacency case
    start_time <- Sys.time()
    ad_vec <- svd(A, nu = K, nv = K);
    ad_vec <- ad_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(ad_vec)%*%Pi%*%M[,j], type = "2"));
    }
    tau <- min(tau_set)l
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "NGfs"] <- tau;
    mj_mat[ii, "NGfs"] <- Mj;
    time_mat[ii, "NGfs"] <- end_time - start_time;
    
    start_time <- Sys.time()
    ad_vec <- svd(A, nu = 2*K, nv = 2*K);
    ad_vec <- ad_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(ad_vec)%*%Pi%*%M[,j], type = "2"));
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
    lp_vec <- svd(Net, nu = K, nv = K);
    lp_vec <- lp_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(lp_vec)%*%Pi%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time()
    
    tau_mat[ii, "NGfs2"] <- tau;
    mj_mat[ii, "NGfs2"] <- Mj;
    time_mat[ii, "NGfs2"] <- end_time - start_time;
    
    start_time <- Sys.time()
    lp_vec <- svd(Net, nu = 2*K, nv = 2*K);
    lp_vec <- lp_vec$u;
    tau_set <- c();
    for (j in signals){
      tau_set<-c(tau_set,norm(t(lp_vec)%*%Pi%*%M[,j], type = "2"));
    }
    tau <- min(tau_set);
    Mj <- norm(M[,signals[which(tau_set == tau)]], type = "2");
    end_time <- Sys.time();
    
    tau_mat[ii, "NGfs22k"] <- tau;
    mj_mat[ii, "NGfs22k"] <- Mj;
    time_mat[ii, "NGfs22k"] <- end_time - start_time;
  }
  
  filename = paste('./Exp1_a_', a, '.Rdata', sep = "");
  save(tau_mat, mj_mat, time_mat, file = filename)
  tau_matrix[ajj,] = colMeans(tau_mat); mj_matrix[ajj,] = apply(mj_mat, 2, mean)
  
  print(paste("a = ", a));
}

save(aseq, tau_matrix, mj_matrix, file = "fig1:left-dcsbm.Rdata")

