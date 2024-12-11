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

cluster <- function(T){
  n = sum(T); k = dim(T)[1]; truek = dim(T)[2];
  if(k > truek) stop("Error: Number of estimated clusters is larger than truth")
  match = 0;
  
  if(sum(T) == 0){
    return(list(error = 0, flagrecord = 1:k))
  }   
  
  if(k == 1) {match = max(T);  flagrecord = which.max(T);}  
  if(k == 2 & truek == 2) {
    m = c(sum(diag(T)),  n - sum(diag(T)) );
    match = max(m); flagrecord = c(which.max(m), 3 - which.max(m));
  }
  
  if(k == 2 & truek == 3){
    for(i1 in 1:3){
      m = 0; flag = numeric(3);
      flag(i1) = 1; 
      m = m + T(1, i1);
      m = m + max(T[2, flag == 0]);
      if(match < m)
      {match = m; flagrecord = c(i1, which(T[2, ] == max(T[2, flag == 0])) );}
    }
  }
  if(k == 3 & truek == 3){
    for(i1 in 1:3){
      m = 0; flag = numeric(3);
      flag[i1] = 1; 
      m = m + T[1, i1];
      for(i2 in 1:3){
        if(flag[i2]==0) 
        {flag[i2] = 1; m = m + T[2, i2];}
        else
          next;
        m = m + T[3, flag == 0];
        if(m > match)
          match = m; flagrecord = c(i1, i2, which(flag==0));
          m = T[1, i1];
          flag[i2] = 0;
      }
    }
  }
  
  if(truek == 4){
    for(i1 in 1:4){
      flag = numeric(4);
      flag[i1] = 1;
      T2 = T[2:k, flag==0];
      result = cluster(T2); 
      m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
      flag[flag >= i1] = flag[flag >= i1] + 1;
      if(match < m) {match = m; flagrecord = c(i1, flag);}
    }
  }
  
  if(truek == 5){
    for(i1 in 1:5){
      flag = numeric(5);
      flag[i1] = 1;
      T2 = T[2:k, flag==0];
      result = cluster(T2); 
      m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
      flag[flag >= i1] = flag[flag >= i1] + 1;
      if(match < m) {match = m; flagrecord = c(i1, flag);}
    }
  }
  
  if(truek == 6){
    for(i1 in 1:6){
      flag = numeric(6);
      flag[i1] = 1;
      T2 = T[2:k, flag==0];
      result = cluster(T2); 
      m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
      flag[flag >= i1] = flag[flag >= i1] + 1;
      if(match < m) {match = m; flagrecord = c(i1, flag);}
    }
  }
  
  if(truek == 7){
    for(i1 in 1:7){
      flag = numeric(7);
      flag[i1] = 1;
      T2 = T[2:k, flag==0];
      result = cluster(T2); 
      m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
      flag[flag >= i1] = flag[flag >= i1] + 1;
      if(match < m) {match = m; flagrecord = c(i1, flag);}
    }
  }
  error = 1 - match/n;    
  return(list(error = error, flagrecord = flagrecord))
}



regp <- function(lm_model){
  f <- summary(lm_model)$fstatistic;
  p <- pf(f[1],f[2],f[3],lower.tail = F);
  return(p)
}

PCAclu <- function(X, K){
  zz = svd(X, nu = K, nv = K);
  Xi = zz$u%*%diag(zz$d[1:K]);
  attempt <- try(expr = {kmeans(Xi, centers = K, nstart = 100)},silent = TRUE)
  estl = kmeans(Xi, K, nstart = 100)
  if(class(attempt) == "try-error"){
    mm = apply(Xi, 1, function(x) sum(x^2));
    nonsparse = which(abs(mm - 1) <= 0.01)
    kcenters = sample(nonsparse, K*2);
    starter = unique(Xi[kcenters,])[1:K,]
    estl = kmeans(Xi, centers = starter)
    best <- sum(estl$wss)
    for(kmeansrep in 2:100){
      kcenters = sample(nonsparse, K*2);
      starter = unique(Xi[kcenters,])[1:K,]
      ee = kmeans(Xi, centers = starter)
      if(sum(ee$wss) < best){
        estl <- ee; best <- sum(ee$wss);
      }
    }
  }
  return(estl$cluster)
}

PCAclu_estK <- function(X, K, estK){
  zz = svd(X, nu = estK, nv = estK);
  Xi = zz$u%*%diag(zz$d[1:estK]);
  attempt <- try(expr = {kmeans(Xi, centers = K, nstart = 100)},silent = TRUE)
  estl = kmeans(Xi, K, nstart = 100)
  if(class(attempt) == "try-error"){
    mm = apply(Xi, 1, function(x) sum(x^2));
    nonsparse = which(abs(mm - 1) <= 0.01)
    kcenters = sample(nonsparse, K*2);
    starter = unique(Xi[kcenters,])[1:K,]
    estl = kmeans(Xi, centers = starter)
    best <- sum(estl$wss)
    for(kmeansrep in 2:100){
      kcenters = sample(nonsparse, K*2);
      starter = unique(Xi[kcenters,])[1:K,]
      ee = kmeans(Xi, centers = starter)
      if(sum(ee$wss) < best){
        estl <- ee; best <- sum(ee$wss);
      }
    }
  }
  return(estl$cluster)
}


sas <- function(X, k,tot, initial_set, s, itermax){
  p = dim(X)[2]
  set0 = initial_set
  set1 = c(0)
  iternum = 0
  rand_index = rep(NA,itermax)
  diff = rep(NA,itermax)
  while(iternum< itermax){
    clustering = kmeans(X[,set0], centers=k,nstart=2)
    result = clustering$cluster
    group = list()
    cond = TRUE
    for(j in seq(1,k)){
      group[[j]] = which(result == j)
      cond = cond && length(group[[j]]) > 1 
    }
    center = NULL
    wcss = rep(0,p)
    if(cond){
      for(j in seq(1,k)){
        center = rbind(center,colMeans(X[group[[j]],]))
        Xc =  t(apply( X[group[[j]],], 1, function(x) x-center[j,]))
        wcss = wcss + apply(Xc, 2, function(x){sum(x^2)})
      }
      iternum = iternum + 1
    }
    set1 = set0
    set0 = which(rank((tot-wcss)/tot,ties.method = "random") > p-s)
    rand_index[iternum] = RRand(true_clust, result)$Rand
    diff[iternum] = length(setdiff(set1,SET)) + length(setdiff(SET,set1))
  }
  out = list(diff = diff, rand_index = rand_index,final_set = set0, iternum = iternum, result = result, betweenss = clustering$betweenss)
  return(out)
}


