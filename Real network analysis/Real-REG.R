source("./Functions/Functions-S1.R")
source("./Real\ network\ analysis/Data.R")
require(logitnorm)
require(ncvreg)

sseq <- seq(50,500,50);
repetition = 50;

## store squared mse
clesina = matrix(0, nrow = length(sseq), ncol = 5);
colnames(clesina) = c("NGfs", "NGfslasso", "lasso",  "mcp", "scad");

clesina_var = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina_var) = c("NGfs", "NGfslasso", "lasso",  "mcp", "scad");

## store the number of positives 
clesina2 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina2) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");

clesina_var2 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina_var2) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");

## store the number of true positives 
clesina3 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina3) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");

clesina_var3 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina_var3) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");


##### S1:sigma = 0.2
for(sjj in 1:length(sseq)){
  s = sseq[sjj];
  
  clesina_mat = matrix(0, nrow = repetition, ncol = 4); 
  colnames(clesina_mat) = c("NGfs", "lasso", "mcp", "scad");
  
  clesina_mat2 = matrix(0, nrow = repetition, ncol = 4); 
  colnames(clesina_mat2) = c("NGfs", "lasso", "mcp", "scad");
  
  clesina_mat3 = matrix(0, nrow = repetition, ncol = 4); 
  colnames(clesina_mat3) = c("NGfs", "lasso", "mcp", "scad");
  
  for (i in 1:repetition){
    net_set <- 1:2000;
    sina_X <- matrix(0, nrow = n, ncol=p);
    for (ii in 1:p){
      s1 = runif(1, -3.5, -3.3);
      h1 = runif(1, 1.5, 2.2);
      sina_X[,ii] = rlogitnorm(n = n, mu = s1, sigma = h1);
    }
    value <- sample(1:p, 10);
    for (ii in 1:10){
      sina_X[,value[ii]] = sina_content[,ii];
    }
    fsresult_K = NGCS(sina_A[net_set,net_set], sina_X[net_set, ], K, Standardize = TRUE, Type = "lsingular"); 
    
    ## Find the response variable we are using, or you can use another one
    vars<-c();
    for (kkk in 1:length(value)) vars <- c(vars, var(sina_X[,value[kkk]]));
    aim_var <- value[which(round(vars,4) == 0.0346)];
    if (length(which(fsresult_K$selected == aim_var)) == 1){
      fs_select <- fsresult_K$selected[-which(fsresult_K$selected==aim_var)];
    }else{
      fs_select <- fsresult_K$selected;
    }
    
    sig_vec = rep(0,p);
    sig_vec2 = rep(0,p-1);
    sig_vec[value] = 1;
    value2 <- c();
    for (ii in value){
      if (ii != aim_var){
        if (ii < aim_var){
          value2 <- c(value2, ii);
        }else{
          value2 <- c(value2, ii-1);
        }
      }
    }
    sig_vec2[value2] = 1;
    
    all_train_set <- (1:n)[-net_set];
    train_set <- sample(all_train_set, s);
    train_X <- sina_X[train_set, aim_var]+rnorm(length(train_set),0,0.2);
    
    data_ngfsx <- as.data.frame(cbind(train_X, sina_X[train_set,fs_select]));
    lm_model <- lm(train_X~., data = data_ngfsx);
    
    lasso_model <- cv.ncvreg(sina_X[train_set,-aim_var], train_X, family = "gaussian", penalty = "lasso");
    
    mcp_model <- cv.ncvreg(sina_X[train_set,-aim_var], train_X, family = "gaussian", penalty = "MCP");
    
    scad_model <- cv.ncvreg(sina_X[train_set,-aim_var], train_X, family = "gaussian", penalty = "SCAD");
    
    
    
    lm_predict <- predict(lm_model, as.data.frame(cbind(train_X,sina_X[net_set, fs_select])));
    lasso_predict <- predict(lasso_model, sina_X[net_set,-aim_var], s = "lambda.min");
    mcp_predict <- predict(mcp_model, sina_X[net_set,-aim_var], s = "lambda.min");
    scad_predict <- predict(scad_model, sina_X[net_set,-aim_var], s = "lambda.min");
    
    clesina_mat[i,"NGfs"] = sqrt(mean((lm_predict - sina_X[net_set,aim_var])^2));
    clesina_mat[i, "lasso"] = sqrt(mean((lasso_predict - sina_X[net_set,aim_var])^2));
    clesina_mat[i, "mcp"] = sqrt(mean((mcp_predict - sina_X[net_set,aim_var])^2));
    clesina_mat[i, "scad"] = sqrt(mean((scad_predict - sina_X[net_set,aim_var])^2));
    
    clesina_mat2[i,"NGfs"] = length(fs_select);
    clesina_mat2[i, "lasso"] = length(which(coef(lasso_model)!=0))-1;
    clesina_mat2[i, "mcp"] = length(which(coef(mcp_model)!=0))-1;
    clesina_mat2[i, "scad"] = length(which(coef(scad_model)!=0))-1;
    
    clesina_mat3[i,"NGfs"] = sum(sig_vec[fs_select]);
    clesina_mat3[i, "lasso"] = sum(sig_vec2[(which(coef(lasso_model)!=0)-1)[-1]]);
    clesina_mat3[i, "mcp"] = sum(sig_vec2[(which(coef(mcp_model)!=0)-1)[-1]]);
    clesina_mat3[i, "scad"] = sum(sig_vec2[(which(coef(scad_model)!=0)-1)[-1]]);
    
    
    
  }
  
  clesina[sjj,] = colMeans(clesina_mat);
  clesina2[sjj,] = colMeans(clesina_mat2);
  clesina3[sjj,] = colMeans(clesina_mat3);
  
  for (t in 1:ncol(clesina)) clesina_var[sjj,t] = sd(clesina_mat[,t]);
  for (t in 1:ncol(clesina2)) clesina_var2[sjj,t] = sd(clesina_mat2[,t]);
  for (t in 1:ncol(clesina3)) clesina_var3[sjj,t] = sd(clesina_mat3[,t]);
  
  print(paste("s = ", s));
}


save(cbind(clesina, clesina2, clesina3, clesina_var, clesina_var2, clesina_var3), file = "Table-1.Rdata")



##### S2:sigma = 0.5

## store squared mse
clesina4 = matrix(0, nrow = length(sseq), ncol = 5);
colnames(clesina4) = c("NGfs", "NGfslasso", "lasso",  "mcp", "scad");

clesina_var4 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina_var4) = c("NGfs", "NGfslasso", "lasso",  "mcp", "scad");

## store the number of positives 
clesina5 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina5) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");

clesina_var5 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina_var5) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");

## store the number of true positives 
clesina6 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina6) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");

clesina_var6 = matrix(0, nrow = length(sseq), ncol = 5); 
colnames(clesina_var6) = c("NGfs", "NGfslasso", "lasso", "mcp", "scad");

for(sjj in 1:length(sseq)){
  s = sseq[sjj];
  
  clesina_mat = matrix(0, nrow = repetition, ncol = 4); 
  colnames(clesina_mat) = c("NGfs", "lasso", "mcp", "scad");
  
  clesina_mat2 = matrix(0, nrow = repetition, ncol = 4); 
  colnames(clesina_mat2) = c("NGfs", "lasso", "mcp", "scad");
  
  clesina_mat3 = matrix(0, nrow = repetition, ncol = 4); 
  colnames(clesina_mat3) = c("NGfs", "lasso", "mcp", "scad");
  
  for (i in 1:repetition){
    net_set <- 1:2000;
    sina_X <- matrix(0, nrow = n, ncol=p);
    for (ii in 1:p){
      s1 = runif(1, -3.5, -3.3);
      h1 = runif(1, 1.5, 2.2);
      sina_X[,ii] = rlogitnorm(n = n, mu = s1, sigma = h1);
    }
    value <- sample(1:p, 10);
    for (ii in 1:10){
      sina_X[,value[ii]] = sina_content[,ii];
    }
    fsresult_K = NGCS(sina_A[net_set,net_set], sina_X[net_set, ], K, Standardize = TRUE, Type = "lsingular"); 
    
    ## Find the response variable we are using, or you can use another one
    vars<-c();
    for (kkk in 1:length(value)) vars <- c(vars, var(sina_X[,value[kkk]]));
    aim_var <- value[which(round(vars,4) == 0.0346)];
    if (length(which(fsresult_K$selected == aim_var)) == 1){
      fs_select <- fsresult_K$selected[-which(fsresult_K$selected==aim_var)];
    }else{
      fs_select <- fsresult_K$selected;
    }
    
    sig_vec = rep(0,p);
    sig_vec2 = rep(0,p-1);
    sig_vec[value] = 1;
    value2 <- c();
    for (ii in value){
      if (ii != aim_var){
        if (ii < aim_var){
          value2 <- c(value2, ii);
        }else{
          value2 <- c(value2, ii-1);
        }
      }
    }
    sig_vec2[value2] = 1;
    
    all_train_set <- (1:n)[-net_set];
    train_set <- sample(all_train_set, s);
    train_X <- sina_X[train_set, aim_var]+rnorm(length(train_set),0,0.5);
    
    data_ngfsx <- as.data.frame(cbind(train_X, sina_X[train_set,fs_select]));
    lm_model <- lm(train_X~., data = data_ngfsx);
    
    lasso_model <- cv.ncvreg(sina_X[train_set,-aim_var], train_X, family = "gaussian", penalty = "lasso");
    
    mcp_model <- cv.ncvreg(sina_X[train_set,-aim_var], train_X, family = "gaussian", penalty = "MCP");
    
    scad_model <- cv.ncvreg(sina_X[train_set,-aim_var], train_X, family = "gaussian", penalty = "SCAD");
    
    
    
    lm_predict <- predict(lm_model, as.data.frame(cbind(train_X,sina_X[net_set, fs_select])));
    lasso_predict <- predict(lasso_model, sina_X[net_set,-aim_var], s = "lambda.min");
    mcp_predict <- predict(mcp_model, sina_X[net_set,-aim_var], s = "lambda.min");
    scad_predict <- predict(scad_model, sina_X[net_set,-aim_var], s = "lambda.min");
    
    clesina_mat[i,"NGfs"] = sqrt(mean((lm_predict - sina_X[net_set,aim_var])^2));
    clesina_mat[i, "lasso"] = sqrt(mean((lasso_predict - sina_X[net_set,aim_var])^2));
    clesina_mat[i, "mcp"] = sqrt(mean((mcp_predict - sina_X[net_set,aim_var])^2));
    clesina_mat[i, "scad"] = sqrt(mean((scad_predict - sina_X[net_set,aim_var])^2));
    
    clesina_mat2[i,"NGfs"] = length(fs_select);
    clesina_mat2[i, "lasso"] = length(which(coef(lasso_model)!=0))-1;
    clesina_mat2[i, "mcp"] = length(which(coef(mcp_model)!=0))-1;
    clesina_mat2[i, "scad"] = length(which(coef(scad_model)!=0))-1;
    
    clesina_mat3[i,"NGfs"] = sum(sig_vec[fs_select]);
    clesina_mat3[i, "lasso"] = sum(sig_vec2[(which(coef(lasso_model)!=0)-1)[-1]]);
    clesina_mat3[i, "mcp"] = sum(sig_vec2[(which(coef(mcp_model)!=0)-1)[-1]]);
    clesina_mat3[i, "scad"] = sum(sig_vec2[(which(coef(scad_model)!=0)-1)[-1]]);
    
    
    
  }
  
  clesina4[sjj,] = colMeans(clesina_mat);
  clesina5[sjj,] = colMeans(clesina_mat2);
  clesina6[sjj,] = colMeans(clesina_mat3);
  
  for (t in 1:ncol(clesina4)) clesina_var4[sjj,t] = sd(clesina_mat[,t]);
  for (t in 1:ncol(clesina5)) clesina_var5[sjj,t] = sd(clesina_mat2[,t]);
  for (t in 1:ncol(clesina6)) clesina_var6[sjj,t] = sd(clesina_mat3[,t]);
  
  print(paste("s = ", s));
}


save(cbind(clesina4, clesina5, clesina6, clesina_var4, clesina_var5, clesina_var6), file = "Table-2.Rdata")
