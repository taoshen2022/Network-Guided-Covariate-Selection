source("./Functions/Functions-S2.R")
require(ncvreg)
require(care)

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

mseq = seq(1.0, 2.0, 0.1)

mse = matrix(0, nrow = length(mseq), ncol = 7) #store the error rate in each repetition
colnames(mse) = c("Lasso", "NGreg", "SLM", "PCR", "MCP", "SCAD", "Base");
mse2 = mse;

#set up the difference between the different communities
for(mujj in 1:length(mseq)){
  mu = mseq[mujj]
  
  ############ End Parameter Set  ###########
  # regression error
  msemat = matrix(0, nrow = repetition, ncol = 7) #store the error rate in each repetition
  colnames(msemat) = c("Lasso", "NGreg", "SLM", "PCR", "MCP", "SCAD", "Base");
  msemat2 = msemat;
  

  for (ii in 1: repetition){
    set.seed(ii+9999)
    
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
    X = Xmean + matrix(rnorm(n*p, 0, sqrt(0.5)), n, p); 
    
    # set up the response variable
    gamma = matrix(abs(rnorm(L, 0, 1)), nrow = L, ncol = 1);
    gamma = gamma - mean(gamma);
    Zmean = Y%*%gamma;
    Z = Zmean + rnorm(n, 0, sqrt(0.5));
    
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
    fsresult = NGCS(A, X[1:n1, ], L, Lap = FALSE, Standardize = FALSE); 
    pval.new = fsresult$pval;
    
    # New approach: Network-Guided, Laplacian
    # fslap = NGCS(A, X[1:n1, ], L, Lap = TRUE, Standardize = FALSE); 
    # pval.lap = fslap$pval;
    
    
    ######### Regression ###################
    # all covariates with Lasso
    cv_model <- cv.glmnet(X[(n1+1):n,], Z[(n1+1):n], alpha = 1);
    best_lambda <- cv_model$lambda.min;
    
    zzfit = glmnet(X[(n1+1):n,], Z[(n1+1):n], alpha = 1, lambda = best_lambda, family = "gaussian");
    coef = coef(zzfit);
    est.lasso = predict(zzfit, newx = X[1:n1,], s = best_lambda);
    msemat[ii, "Lasso"] = mean((est.lasso - Z[1:n1])^2);
    msemat2[ii, "Lasso"] = mean((est.lasso - Zmean[1:n1])^2);
    
    # Regression with NGFS
    est.fs = reg_pred(X[(n1+1):n,fsresult$selected], X[1:n1,fsresult$selected], Z[(n1+1):n], min(L, fsresult$S));
    msemat[ii, "NGreg"] = mean((est.fs - Z[1:n1])^2);
    msemat2[ii, "NGreg"] = mean((est.fs - Zmean[1:n1])^2);
    
    
    # Regression with NGFS, Laplacian
    # est.fslap = PCAclu(X[,fslap$selected], K);
    # errormat[ii, "NG clustering 2"] = cluster(table(est.fslap, l))$error;
    # nmimat[ii, "NG clustering 2"] = NMI(est.fslap, l);
    
    # PCR
    est.pcr <- reg_pred(X[(n1+1):n,], X[1:n1,], Z[(n1+1):n], L);
    msemat[ii, "PCR"] = mean((est.pcr - Z[1:n1])^2);
    msemat2[ii, "PCR"] = mean((est.pcr - Zmean[1:n1])^2);
    
    
    # CAR
    car <- carscore(X[(n1+1):n,], Z[(n1+1):n], lambda=0);
    car.pval = 1-pbeta(car^2, shape1=1/2, shape2=(n1-2)/2);
    car.re = HCfs(car.pval);
    car.fit = slm(X[(n1+1):n,car.re$selected], Z[(n1+1):n], lambda=0, lambda.var=0);
    est.car <- predict(car.fit, as.matrix(X[1:n1,car.re$selected]));
    msemat[ii, "SLM"] = mean((est.car - Z[1:n1])^2);
    msemat2[ii, "SLM"] = mean((est.car - Zmean[1:n1])^2, na.rm = TRUE);
    
    # MCP SCAD
    mcp_model <- cv.ncvreg(X[(n1+1):n,], Z[(n1+1):n], family = "gaussian", penalty = "MCP");
    scad_model <- cv.ncvreg(X[(n1+1):n,], Z[(n1+1):n], family = "gaussian", penalty = "SCAD");
    mcp_predict <- predict(mcp_model, X[1:n1,], s = "lambda.min");
    scad_predict <- predict(scad_model, X[1:n1,], s = "lambda.min");
    
    msemat[ii, "MCP"] = mean((mcp_predict - Z[1:n1])^2);
    msemat2[ii, "MCP"] = mean((mcp_predict - Zmean[1:n1])^2, na.rm = TRUE);
    msemat[ii, "SCAD"] = mean((scad_predict - Z[1:n1])^2);
    msemat2[ii, "SCAD"] = mean((scad_predict - Zmean[1:n1])^2, na.rm = TRUE);
    
    
    msemat[ii, "Base"] = mean((mean(Z[1:n1]) - Z[1:n1])^2);
    msemat2[ii, "Base"] = mean((mean(Z[1:n1]) - Z[1:n1])^2);
    
    print(ii)
  }
  
  mse[mujj,] = colMeans(msemat); mse2[mujj,] = colMeans(msemat2); 
  
  print(paste("mu = ", mu));
}


save(mseq, mse, file = "fig2-right.Rdata")
