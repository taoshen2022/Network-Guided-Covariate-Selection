require(RSpectra)
require(igraph)
require(Rdimtools)
require(factoextra)
require(cluster)
require(pracma)
require(mclust)
require(phyclust)
require(Matrix)
require(ggplot2)
require(grid)
require(gtable)
require(mvtnorm)

NGCS <- function(A, X, K, Lap = TRUE, Standardize = TRUE, S = NULL, alpha = 1/2, pvalcut = NULL, Type = "eigen"){
  #check the dimensions of A and X
  n = dim(X)[1]; p = dim(X)[2]
  if(dim(A)[1]!=n) stop("Error: the dimension of A does not match the dimension of X!")
  if(dim(A)[2]!=n) stop("Error: A should be a square matrix!")
  
  Net = A;
  if(Lap){#define the Laplacian matrix for the network information
    d = colSums(A); #degree vector
    ind_neighbor = which(d > 0); #indices for nodes with at least one 
    Dtau=d[ind_neighbor]^(-0.5);
    Ltau = t(Dtau * A[ind_neighbor, ind_neighbor]) * Dtau;
    Net[ind_neighbor, ind_neighbor] = Ltau;
  }
  
  if(Standardize){#standardize the covariates so that each with zero mean and standard deviation 1
    X = apply(X, 2, function(x) x - mean(x));
    X = apply(X, 2, function(x) sqrt(n-1)*x/sqrt(sum(x^2)));
  }
  
  if (Type == "eigen"){
    Net_eigen = RSpectra::eigs(Net, k = K);
    vec = Net_eigen$vectors; 
    
    t = (vec[,1]%*%X)^2;
    if(K >= 2){
      for(kk in 2:K){
        t = t + (vec[,kk]%*%X)^2;
      } 
    }
    #t = as.numeric(t)
    pval = pchisq(t, df = K, lower.tail = FALSE);
    if(!is.null(S)){
      psort = sort(pval, decreasing = FALSE);
      selected = which(pval <= psort[S]) #It will keep all the ties
      return(list(pval = pval, selected = selected, hc = NULL, S = S, K = K))
    }
  }
  if (Type == "lsingular"){
    Net_svd = RSpectra::svds(Net, K, nu = K, nv = K)
    Net_u = Net_svd$u
    Net_v = Net_svd$v
    #vec = Net_eigen$vectors; 
    
    t = (Net_u[,1]%*%X)^2;
    if(K >= 2){
      for(kk in 2:K){
        t = t + (Net_u[,kk]%*%X)^2;
      } 
    }
    #t = (Net_v[,1]%*%X)^2;
    #if(K >= 2){
    #for(kk in 2:K){
    #t = t + (Net_v[,kk]%*%X)^2;
    #} 
    #}
    #t = as.numeric(t)
    pval = pchisq(t, df = K, lower.tail = FALSE);
    if(!is.null(S)){
      psort = sort(pval, decreasing = FALSE);
      selected = which(pval <= psort[S]) #It will keep all the ties
      return(list(pval = pval, selected = selected, hc = NULL, S = S, K = K))
    }
  }
  if (Type == "rsingular"){
    
    Net_svd = RSpectra::svds(Net, K, nu = K, nv = K)
    Net_u = Net_svd$u
    Net_v = Net_svd$v
    #vec = Net_eigen$vectors; 
    
    #t = (Net_u[,1]%*%X)^2;
    #if(K >= 2){
    #for(kk in 2:K){
    #t = t + (Net_u[,kk]%*%X)^2;
    #} 
    #}
    t = (Net_v[,1]%*%X)^2;
    if(K >= 2){
      for(kk in 2:K){
        t = t + (Net_v[,kk]%*%X)^2;
      } 
    }
    #t = as.numeric(t)
    pval = pchisq(t, df = K, lower.tail = FALSE);
    if(!is.null(S)){
      psort = sort(pval, decreasing = FALSE);
      selected = which(pval <= psort[S]) #It will keep all the ties
      return(list(pval = pval, selected = selected, hc = NULL, S = S, K = K))
    }
  }
  if (Type == "singular"){
    
    Net_svd = RSpectra::svds(Net, K, nu = K, nv = K)
    Net_u = Net_svd$u
    Net_v = Net_svd$v
    #vec = Net_eigen$vectors; 
    
    t = (Net_u[,1]%*%X)^2;
    if(K >= 2){
      for(kk in 2:K){
        t = t + (Net_u[,kk]%*%X)^2;
      } 
    }
    t = t + (Net_v[,1]%*%X)^2;
    if (K >= 2){
      for(kk in 2:K){
        t = t + (Net_v[,kk]%*%X)^2;
      } 
    }
    #t = as.numeric(t)
    pval = pchisq(t, df = K, lower.tail = FALSE);
    if(!is.null(S)){
      psort = sort(pval, decreasing = FALSE);
      selected = which(pval <= psort[S]) #It will keep all the ties
      return(list(pval = pval, selected = selected, hc = NULL, S = S, K = K))
    }
  }
  
  
  HCresult = HCfs(pval, alpha = alpha, pvalcut = pvalcut);
  HC = HCresult$HC;
  S = HCresult$S; 
  selected = HCresult$selected
  return(list(pval = pval, selected = selected, hc = HC, S = S, K = K))
  
}


HCfs <- function(pval, alpha = 1/2, pvalcut = NULL){
  # with given p-values, decide whether there is a signal or not
  
  n = length(pval);
  if(is.null(pvalcut)){pvalcut = 2/n;}
  
  kk = (1:n)/(1 + n);
  psort = sort(pval, decreasing = FALSE);
  ind = which(psort > pvalcut & kk <= alpha);
  HC = sqrt(n)*(kk[ind] - psort[ind])/sqrt(psort[ind] - psort[ind]^2);
  HCT = max(HC);
  z <- which(HC == HCT); 
  numselect = ind[z[length(z)]];
  select = which(pval <= psort[numselect])
  
  return(list(selected = select, HC = HC, S = numselect));
}

HCfs_cri <- function(pval, alpha = 1/2, cri){
  # with given number cri, we select cri p-values to decide
  
  n = length(pval);
  
  kk = (1:n)/(1 + n);
  psort = rank(pval);
  select = which(psort <= cri)
  HC = sqrt(n)*(kk[select] - pval[select])/sqrt(pval[select] - pval[select]^2);
  
  return(list(selected = select, HC = HC));
}


HCfs_cluster <- function(pval, samplesize, alpha = 1/2, pvalcut = NULL){
  # with given p-values, decide whether there is a signal or not
  
  n = length(pval);
  if(is.null(pvalcut)){pvalcut = 2/n;}
  
  kk = (1:n)/(1 + n);
  psort = sort(pval, decreasing = FALSE);
  ind = which(psort > pvalcut & kk <= alpha);
  HC = sqrt(n)*(kk[ind] - psort[ind])/sqrt(kk[ind]);
  HC = HC/sqrt(max(sqrt(samplesize)*(kk[ind] - psort[ind])/kk[ind], 0) + 1 );
  HCT = max(HC);
  z <- which(HC == HCT); 
  numselect = ind[z[length(z)]];
  select = which(pval <= psort[numselect])
  
  
  return(list(selected = select, HC = HC, S = numselect));
}






regp <- function(lm_model){
  f <- summary(lm_model)$fstatistic;
  p <- pf(f[1],f[2],f[3],lower.tail = F);
  return(p)
}

reg_pred <- function(X, newX, z, K){
  Xsvd = svd(X, nu = K, nv = K);
  z = matrix(z, nrow = 1, ncol = dim(X)[1])
  coef = z%*%(Xsvd$u%*%diag(1/Xsvd$d[1:K])%*%t(Xsvd$v));
  pred = newX%*%t(coef);
  return(pred)
}


latent.vector<-function(n, L, mu = 0){
  num_m <- c()
  p = runif(n*L)+1
  p2 = runif(n*L)+1
  for (i in 1:(n*L)){
    if (p[i] <= 1/2){
      num_m <- c(num_m, p2[i]^2)
    } else {
      num_m <- c(num_m, (1-p2[i])^2)
    }
  } 
  m = matrix(num_m, nrow = n, ncol = L);
  return(m)
}






