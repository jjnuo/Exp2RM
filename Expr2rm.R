rm(list=ls())
setwd('/XXXXXXXXXXXXXXXXXXXXXXXXXXXX/')
options(digits=5) 
library(glmnet) 
library(readr)
library(dplyr)
library(doMC)
registerDoMC(cores = 4) #parallel computing with 4 cores
set.seed(123)

rsq <- function(pred,actual){
  # Sum of squares total and error
  sst <- sum((actual - mean(actual))^2)
  sse <- sum((actual - pred)^2)
  # R squared
  rsq <- 1 - sse / sst
}

train <- function (){
  for (n in seq_len(site_number)){
    set.seed(123)
    best_lambda<-list()
    best_mse <-list()
    best_rs<-list()
    min_or_1se<-list()
    for (i in 1:9){
      set.seed(123)
      alpha= i/10
      cv.fit<-cv.glmnet(Xtrain, Ytrain[,n,drop = FALSE],
                        nfolds = nfold,
                        n_lambda=n_lambda,
                        grouped=FALSE, 
                        alpha= alpha, 
                        type.measure = "mse",
                        parallel = TRUE, #parallel computing
                        keep = TRUE)
      
      
      # Compare lambda.min and lambda.1se
      ass.min<- assess.glmnet(cv.fit, 
                              newx = Xtest, newy = Ytest[,n,drop=FALSE], 
                              s="lambda.min",
                              family="gaussian")
      ass.1se<- assess.glmnet(cv.fit, 
                              newx = Xtest, newy = Ytest[,n,drop=FALSE], 
                              s="lambda.1se",
                              family="gaussian")
      
      mse.min <- ass.min[[1]][1] #mse of the predicts
      mse.1se <- ass.1se[[1]][1] #mse of the predicts
      # choose either lambda.min or lambda.1se based on the mse of predicts
      if (mse.min <= mse.1se ){
        best_mse<-append(best_mse,cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.min])
        best_lambda<-append(best_lambda, cv.fit$lambda.min)
        best_rs <- append(best_rs,cv.fit$glmnet.fit$dev.ratio[which(cv.fit$glmnet.fit$lambda == cv.fit$lambda.min)])
        min_or_1se<- append(min_or_1se, 1)
        
      } else {
        best_mse<-append(best_mse,cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.1se])
        best_lambda<-append(best_lambda, cv.fit$lambda.1se)
        best_rs <- append(best_rs,cv.fit$glmnet.fit$dev.ratio[which(cv.fit$glmnet.fit$lambda == cv.fit$lambda.1se)])#dev.ratio=1-dev/nulldev.
        min_or_1se<- append(min_or_1se, 2)
      }
      
      message("alpha ",i," done.")
      
    }
    bt<-which.min(unlist(best_mse)) 
    tune[n,1]<-(bt)/10
    tune[n,2]<-best_lambda[[bt]]
    tune[n,3]<-min_or_1se[[bt]]
    tune[n,4]<-best_mse[bt][[1]]
    tune[n,5]<-best_rs[[bt]]
    
    best_model<-glmnet(Xtrain, Ytrain[,n,drop = FALSE],
                       nfolds = nfold,
                       lambda = best_lambda[[bt]],
                       grouped=FALSE, 
                       alpha= (bt)/10, 
                       type.measure = "mse",
                       parallel = TRUE, #parallel computing
                       keep = TRUE)
    
    pred <- predict(best_model , s = best_lambda[[bt]], newx = Xtest)
    ass<- assess.glmnet(best_model, 
                        newx = Xtest, newy = Ytest[,n,drop=FALSE], 
                        s=best_lambda[[bt]],
                        family="gaussian")
    
    tune[n,6] <- ass[[1]][1] #mse of the predicts
    Ypred[,n]<- as.matrix(pred)
    test_rsq<-rsq(pred,Ytest[,n,drop=FALSE]) #rsq of the predicts
    tune[n,7] <- test_rsq
    
    #Gene coefficient
    el.coef<-coef(best_model)
    coef.df<-data.frame(GeneName = el.coef@Dimnames[[1]][el.coef@i + 1][-1], Coef = el.coef@x[-1]) #exclude intercept
    coef.df[,'Label']<-paste0('_site_',n)
    gene_coef<- rbind(gene_coef, coef.df)
    message("_________site ",n," done.")
  }
  
  tissue_rsq[1,1]<-rsq(Ypred[1,,drop=FALSE],Ytest[1,,drop=FALSE])
  tissue_rsq[2,1]<-rsq(Ypred[2,,drop=FALSE],Ytest[2,,drop=FALSE])
  saveRDS(tune,paste0('best_tune','.rds'))
  saveRDS(gene_coef, paste0('gene_coef','.rds'))
  saveRDS(tissue_rsq,paste0('tissue_rsq','.rds'))
  saveRDS(Ypred,paste0('Ypred','.rds'))
}


  Xtrain<-readRDS()
  Xtest<-readRDS()
  
  Ytrain<-readRDS()
  Ytest<-readRDS()
  
  n_lambda = 200
  site_number = ncol(Ytrain)
  nfold = 10
  gene_coef<-list()
  tune<-matrix(as.numeric(NA), ncol = 7, nrow = ncol(Ytrain))
  colnames(tune)<-c('alpha','lambda','min_or_1se','mse','train_rsq','test_mse','test_rsq')  
  gene_coef<-matrix(NA, ncol = 3, nrow = 1)
  colnames(gene_coef)<-c('GeneName','Coef','Label') 
  tissue_rsq<-matrix(NA,nrow=2,ncol=1)
  rownames(tissue_rsq)<-rownames(Ytest)
  colnames(tissue_rsq)<'tissue_rsq'
  Ypred<-matrix(NA,nrow=2,ncol=ncol(Ytest))
  rownames(Ypred)<-rownames(Ytest)
  
  
  Xtrain <- as.matrix(Xtrain)
  Ytrain <- as.matrix(Ytrain)
  Xtest <- as.matrix(Xtest)
  Ytest <- as.matrix(Ytest)
  
  train_model<-train()

















