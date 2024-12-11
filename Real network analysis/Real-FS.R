source("./Functions/Functions-S1.R")
source("./Real\ network\ analysis/Data.R")
require(logitnorm)
require(sparcl)

####Task 1: use NGCS to select variables
#### target plot: estK vs recovered true positives
kseq <- seq(1,50,1);
repetition <- 50;
n <- nrow(sina_content)
n1 <- 2000;
TPsina = matrix(0, nrow = length(kseq), ncol = 2); #store the true positives
colnames(TPsina) = c("TP", "P");

for(kjj in 1:length(kseq)){
  k = kseq[kjj];
  ############ End Parameter Set  ############
  # feature selection error
  TPmat = matrix(0, nrow = repetition, ncol = 2); #store the #TP in each repetition
  colnames(TPmat) = c("TP", "P");
  
  #p = 1000; 
  node_selection <- matrix(0, nrow = repetition, ncol = n1);
  for (ii in 1: repetition){
    node_selection[ii,] <- sort(sample(1:n, n1));
  }
  
  for (ii in 1: repetition){
    ######### Covariate selection ###################
    #Truth
    sig_vec = matrix(0, nrow = 1, ncol = p);
    sig_vec = c(sig_vec);
    
    sina_X <- matrix(0, nrow = n, ncol=p);
    value <- sample(1:p, 10);
    for (i in 1:p){
      s = runif(1, -3.5, -3.3);
      h = runif(1, 1.5, 2.2);
      sina_X[,i] = rlogitnorm(n = n, mu = s, sigma = h);
    }
    for (i in 1:10){
      sina_X[,value[i]] = sina_content[,i];
    }
    sig_vec[value] = 1;
    
    # New approach: Network-Guided
    fsresult = NGCS(sina_A[c(node_selection[ii,]),c(node_selection[ii,])], sina_X[c(node_selection[ii,]),], k,  Standardize = TRUE, Type = "lsingular"); 
    #pval_fs = fsresult$pval;
    TPmat[ii, "TP"] = sum(sig_vec[fsresult$selected]);
    TPmat[ii, "P"] = length(fsresult$selected)
    
  }
  
  TPsina[kjj,] = colMeans(TPmat);
  
  print(paste("estK = ", k));
}


####Task 2: use NGCS to select variables
#### target plot: positives vs recovered true positives

K = 20;

pseq = seq(1,40, by = 2); 

TPsina2 = matrix(0, nrow = length(pseq), ncol = 4); #store the true positives
colnames(TPsina2) = c("NGfs", "KS", "SFS", "Oracle");

pval_table_fs <- matrix(0, nrow = repetition, ncol = p);
pval_table_chi <- matrix(0, nrow = repetition, ncol = p); 
pval_table_sfs <- matrix(0, nrow = repetition, ncol = p); 
pval_table_or <- matrix(0, nrow = repetition, ncol = p); 

for (i in 1:repetition){
  node_select <- sort(sample(1:n, n1));
  
  sina_X <- matrix(0, nrow = n, ncol=p);
  value <- sample(1:p, 10);
  value_table[,i] <- value;
  for (j in 1:p){
    s = runif(1, -3.5, -3.3);
    h = runif(1, 1.5, 2.2);
    sina_X[,j] = rlogitnorm(n = n, mu = s, sigma = h);
  }
  for (j in 1:10){
    sina_X[,value[j]] = sina_content[,j];
  }
  
  ##NGFS
  fsresult_K = NGCS(sina_A[node_select,node_select], sina_X[node_select, ], K, Standardize = TRUE, Type = "lsingular"); 
  pval_table_fs[i,] = fsresult_K$pval;
  
  ###KS
  rep = 100 * p;
  KSvalue = rep(0, rep); kk = (0:n)/n;
  for (ii in 1:rep){
    
    x = rnorm(n); 
    z = (x - mean(x))/sd(x);
    z = z/sqrt(1 - 1/n);
    pi = pnorm(z);
    pi = sort(pi);
    
    KSvalue[ii] = max(max(abs(kk[1:n] - pi)), max(abs(kk[2:(n+1)] - pi)));
  }
  KSvalue = KSvalue*sqrt(n);
  KSvalue = sort(KSvalue);
  KSmean = mean(KSvalue);
  KSstd = sd(KSvalue);
  per = 0.1
  KSmean = mean(KSvalue[1:round(rep*per)]);
  KSstd = sd(KSvalue[1:round(rep*per)]);
  
  kk = (0:n)/n;
  KS = rep(0,p)
  for (j in 1:p){
    pi = pnorm(sina_X[,j]/sqrt(1 - 1/n));
    pi = sort(pi);
    KS[j] = sqrt(n)*max(max(abs(kk[1:n] - pi)), max(abs(kk[2:(n+1)] - pi)));
  }
  
  KS = sort(KS);
  KSm = mean(KS); KSs = sd(KS);
  KSm = mean(KS[1:round(per*p)]); KSs = sd(KS[1:round(per*p)]);
  KS = (KS - KSm)/KSs*KSstd + KSmean;
  
  if_pval = rep(0,p);
  for (ii in 1:p){
    if_pval[ii] = mean(KSvalue > KS[ii]);
  }
  
  pval_table_chi[i,] = if_pval;
  
  # SFS
  skmeans = KMeansSparseCluster(sina_X[node_select,], K)
  pval_table_sfs[i,] = skmeans[[length(skmeans)]]$ws
  
  # Oracle
  pval.labe = apply(sina_X[node_select, ], 2, function(x)  regp(lm(x~as.factor(sina_true[node_select])-1)));
  pval_table_or[i,] = pval.labe
  
}


for(pjj in 1:length(pseq)){
  
  cri = pseq[pjj];
  
  TPmat2 = matrix(0, nrow = reptition, ncol = 4); #store the #TP in each repetition
  colnames(TPmat2) = c("NGfs", "KS", "SFS", "Oracle");
  
  
  #Calculation of the error 
  for (ii in 1:repetition){
    
    #Truth
    sig_vec = matrix(0, nrow = 1, ncol = p)
    sig_vec = c(sig_vec)
    value <- value_table[,ii]
    sig_vec[value] = 1
    
    # NGCS
    fsresult = HCfs_cri(pval_table_fs[ii,], cri = cri)
    TPmat2[ii, "NGfs"] = sum(sig_vec[fsresult$selected]);
    
    # KS
    chiresult = HCfs_cri(pval_table_chi[ii,], cri = cri);
    TPmat2[ii, "KS"] = sum(sig_vec[chiresult$selected]);
    
    # SFS
    skresult <- which(rank(pval_table_sfs[ii,])>=p-cri+1)
    TPmat2[ii, "SFS"] = sum(sig_vec[skresult]);
    
    # Oracle
    oracle = HCfs_cri(pval_table_or[ii,], cri = cri);
    TPmat2[ii, "Oracle"] = sum(sig_vec[oracle$selected]);
    
  }
  
  TPsina2[pjj,] = colMeans(TPmat2);
  
  print(paste("p = ", cri));
}

save(kseq, TPsina, file = "fig3-left.Rdata")
pseq_sina <- pseq
save(pseq_sina, TPsina2, file = "fig3-right.Rdata")




