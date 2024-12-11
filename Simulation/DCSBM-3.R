source("./Functions/Functions-S1.R")
require(sparcl)

n = 1000; 
n1 = 800;
n2 = 200;
repetition = 50;

K = 3; #number of communities
signal_size = 50;

ppseq = as.integer(exp(seq(6,12,0.5)));
error = matrix(0, nrow = length(ppseq), ncol = 6); #store the error rate in each repetition
colnames(error) = c("All covariates", "NG clustering", "NG clustering 2", "IF-PCA", "SKmeans", "SAS");
nmi1 = error; #store the nmi rate in each repetition
diserror = error; disnmi = nmi1;



#set up the difference between the different communities
for(pjj in 1:length(ppseq)){
  p = ppseq[pjj];
  signals = sample(p, signal_size);
  mu = 0.3;
  de = 0.2;
  
  
  ############ End Parameter Set  ############
  errormat = matrix(0, nrow = repetition, ncol = 6); #store the error rate in each repetition
  colnames(errormat) = c("All covariates", "NG clustering", "NG clustering 2", "IF-PCA", "SKmeans", "SAS");
  nmimat = errormat; #store the nmi rate in each repetition
  
  
  #Calculation of the error 
  for (ii in 1: repetition){
    set.seed(ii+9999)
    
    M = matrix(0, nrow = K, ncol = p);
    for (i in signals){
      M[,i] = matrix(rnorm(K, mu, de), nrow = K, ncol = 1);
    }
    M = M*matrix(rbinom(K*p, 1, 1/2)*2-1, nrow = K, ncol = p);
    
    a = pi/4;
    B = matrix(c(sin(a)^2, cos(a)^2/2, cos(a)^2/2, cos(a)^2/2, sin(a)^2, cos(a)^2/2, cos(a)^2/2, cos(a)^2/2, sin(a)^2), ncol=3, byrow = TRUE);
    theta = runif(n); 

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
    A = A[1:n1, 1:n1]; #only the first n1 nodes are observed with the network
    
    
    ######### Covariate selection ###################
    #Truth
    sig_vec = rep(0, p); sig_vec[signals] = 1;
    
    # New approach: Network-Guided
    fsresult = NGCS(A, X[1:n1, ], K, Lap = FALSE, Standardize = FALSE); 
    pval.new = fsresult$pval;
    
    # New approach: Network-Guided, Laplacian
    fslap = NGCS(A, X[1:n1, ], K, Lap = TRUE, Standardize = FALSE); 
    pval.lap = fslap$pval;
    
    # Chisquare
    pval.chi = apply(X, 2, function(x) pchisq(sum(x^2), n, lower.tail = FALSE));
    chiresult = HCfs(pval.chi);
    
    # Oracle
    pval.label = apply(X[1:n1, ], 2, function(x)  regp(lm(x~as.factor(l[1:n1])-1)));
    oracle = HCfs(pval.label);
    
    ######### Clustering  ###################
    # all covariates without selection
    est.all = PCAclu(X, K);
    errormat[ii, "All covariates"] = cluster(table(est.all, l))$error;
    nmimat[ii, "All covariates"] = NMI(est.all, l);
    
    
    # Clustering with NGFS
    est.fs = PCAclu(X[,fsresult$selected], K);
    errormat[ii, "NG clustering"] = cluster(table(est.fs, l))$error;
    nmimat[ii, "NG clustering"] = NMI(est.fs, l);
    
    
    # Clustering with NGFS, Laplacian
    est.fslap = PCAclu(X[,fslap$selected], K);
    errormat[ii, "NG clustering 2"] = cluster(table(est.fslap, l))$error;
    nmimat[ii, "NG clustering 2"] = NMI(est.fslap, l);
    
    # Clustering with IF-PCA
    chiresult = HCfs_cluster(pval.chi, n)
    est.chi = PCAclu(X[,chiresult$selected], K);
    errormat[ii, "IF-PCA"] = cluster(table(est.chi, l))$error;
    nmimat[ii, "IF-PCA"] = NMI(est.chi, l);
    
    
    #sparse k-means clustering
    skmeans = KMeansSparseCluster(X, K)
    skmeans.re = skmeans[[length(skmeans)]]$Cs
    errormat[ii, "SKmeans"] = cluster(table(skmeans.re, l))$error;
    nmimat[ii, "SKmeans"] = NMI(skmeans.re, l);
    
    #sas
    center0 = colMeans(X)
    true_clust = l
    SET = seq(1,50)
    s = 50
    Xc0 = t(apply(X,1,function(x) x-center0))
    tot = apply(Xc0, 2, function(x) {sum(x^2)})
    wcss = rep(0,p)
    #Initialize the important feature set.
    for(j in 1:p){
      clustering = kmeans(X[,j], centers = K)
      wcss[j] = clustering$tot.withinss
    }
    rank0 = rank((tot-wcss)/tot,ties.method = "random")
    initial_set = which(rank0 > p-s)
    out = sas(X, k=3,tot, initial_set, s, itermax=10)
    sas.re <- out$result
    errormat[ii, "SAS"] = cluster(table(sas.re, l))$error;
    nmimat[ii, "SAS"] = NMI(sas.re, l);
    
    print(ii)
  }
  
  
  
  filename = paste('./Exp4_cluster_', p, '.Rdata', sep = "");
  save(errormat, nmimat,file = filename)
  error[pjj,] = colMeans(errormat); nmi1[pjj,] = apply(nmimat, 2, mean)
  
  print(paste("p = ", p));
}


save(pseqq, error, file = "fig2-left.Rdata")
