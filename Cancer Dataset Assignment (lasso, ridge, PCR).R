## assignment statistical methods for bioinformatics
#setwd("C:/Users/Administrator/Documents/OneDrive/R_MyFiles/")
load("VIJVER.Rdata")

#install.packages(c("glmnet", "pls","ROCR")) 
library(ROCR)
library(glmnet)
library(pls)
library(GGally)
library(ggplot2)

###########################################################################
# PART 1 : EXPLORE THE DATA
###########################################################################


summary(data)
dim(data)
ncol(data)


#--- pick a few genes at random and evaluate the association with the response
set.seed(1)


rndCols<-sample(seq(2:ncol(data)), 5, replace = FALSE)
rndTestData<-cbind(meta=as.factor(data[,"meta"]),data[,rndCols])
ggpairs(rndTestData)

lr_mod<-glm(meta~.,data = rndTestData,family = 'binomial')
summary(lr_mod)
anova(lr_mod, test="Chisq")


#--- check for collinearity
corMat<-cor(data[,-1])
(sum(abs(corMat) > 0.8)-4948)/2  #don't count diagonal, and mirror correlation


which(corMat > 0.8 & corMat <0.99, arr.ind = T)
plot(data[,2015],data[,6],xlab = names(data)[2015], ylab = names(data)[6])

###########################################################################
# PART B : CONSTRUCT PREDICTORS
###########################################################################

# 0. preparation : split train <-> test ...
#--------------------------------------------------------------------------
data$meta <- ifelse(data$meta =="DM",1,0)  #DM = 1
set.seed (1)
train = sample (c(TRUE , FALSE), nrow(data), rep =TRUE,prob=c(0.5,0.5)) #
test =(! train )
sum(train)


# 1. LASSO 
#--------------------------------------------------------------------------

#--- 1.1 LASSO fit on train
x= model.matrix (meta~., data )[,-1]
y= data$meta
set.seed (1) 
lasso.mod =cv.glmnet(x[train,], y[ train ], alpha =1,family="binomial")
plot(lasso.mod);title("Cross Validation lasso")
bestlam <- lasso.mod$lambda.min
idxMinLambda <-which(lasso.mod$lambda == bestlam)
lasso.mod$nzero[idxMinLambda]  
lasso.mod$cvm[idxMinLambda] # mean cross validation error 1.230502

#--- 1.2 LASSO test performance 
lasso.pred=predict(lasso.mod,s=bestlam,newx=x[test,],type="response")
lasso.pred <- ifelse(lasso.pred > 0.5,1,0)
lasso.misClasificError <- mean(lasso.pred != y[test])  #37%


#--- 1.3 LASSO fit on all data with best lambda and get coefficients
lasso.mod.full <- glmnet(x,y, alpha =1,family="binomial") #fit with all the data (don't use bestlambda here ! see docs)
lasso.coef= predict (lasso.mod.full, s=bestlam, type ="coefficients")
plot(lasso.coef[-1]);title("coefficients Lasso fit",cex = 0.5)
length(lasso.coef[lasso.coef!=0])  # The  nonzero coefficients  (including intercept)
names(data[which(lasso.coef!=0)])[-1]  #get genes, leave out intercept
lasso.coef[which(lasso.coef!=0)][-1] 
data.frame(genes= names(data[which(lasso.coef!=0)])[-1], coefficients=lasso.coef[which(lasso.coef!=0)][-1])
plot(lasso.mod.full)
      
#--- 1.4 LASSO plot performance ROC
pr <- prediction(lasso.pred, y[test])
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf);title("ROC curve Lasso Fit")  #plot ROC curve
auc.lasso <- performance(pr, measure = "auc")
auc.lasso <- auc.lasso@y.values[[1]]  # area under curve of
auc.lasso

#2 RIDGE
#----------------------------------------------------------------------------
#--- 2.1 RIDGE fit on train
set.seed (1) 
ridge.mod =cv.glmnet(x[train,], y[ train ], alpha =0,family="binomial")
plot(ridge.mod)
bestlam <- ridge.mod$lambda.min
idxMinLambda <-which(ridge.mod$lambda == bestlam)
ridge.mod$nzero[idxMinLambda]   #4948 non zero coefficients
ridge.mod$cvm[idxMinLambda] # mean cross validation error 1.30119

#--- 2.2 RIDGE test performance 
ridge.pred=predict(ridge.mod,s=bestlam,newx=x[test,],type="response")
ridge.pred <- ifelse(ridge.pred > 0.5,1,0)
ridge.misClasificError <- mean(ridge.pred != y[test]) 


#--- 2.3 RIDGE fit on all data with best lambda and get coefficients
ridge.mod.full <- glmnet(x,y, alpha =0,family="binomial") #fit with all the data (don't use bestlambda here ! see docs)
ridge.coef= predict(ridge.mod.full, s=bestlam, type ="coefficients")
plot(ridge.coef[-1]);title("coefficients Ridge fit",cex = 0.5)
ridge.coef.sort <- data.frame(genes = names(data[,-1]),coefficients=abs(ridge.coef[-1]))
sum(ridge.coef!=0)



#--- 2.4 RIDGE plot performance ROC
pr <- prediction(ridge.pred, y[test])
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf);title("ROC curve Ridge Fit")  #plot ROC curve
auc.ridge <- performance(pr, measure = "auc")
auc.ridge <- auc.ridge@y.values[[1]]  # area under curve 
auc.ridge



# 3. PCR : 
#--------------------------------------------------------------------------
#--- 3.1  preparation data, folds

# problemo,werkt niet voor bionomial respons (?)
# cv.pcr =pcr(meta~.,data=data,family=binomial(link=logit),validation="CV")

PcaX <-prcomp(data[train,-1])$x #every column is the PC
DataPCR <- data.frame(cbind("meta"=data[train,"meta"],PcaX) ) #the transformed traindata

k =10 #10 fold cv
nrow = nrow(DataPCR)
set.seed (1)
folds = sample (1:k,nrow,replace = TRUE)

#--- 3.2  fit model 10 folds and for every PCA

cv.errors =matrix(NA,k,nrow, dimnames = list(NULL, paste (1:nrow) ))
for (pcIter in 1:nrow) {
  for (fold in 1:k){
  f <- as.formula(paste("meta~",paste('PC', 1:pcIter,sep="",collapse = "+")))
  lr.mod.pca <- glm(f,data = DataPCR[folds!=fold,],family = 'binomial')
  pred= predict (lr.mod.pca, DataPCR[folds ==fold,], id=pcIter,type="response")
  pred<-ifelse(pred > 0.5 ,1,0)
  cv.errors [fold,pcIter]= mean( (DataPCR$meta[folds ==fold]!= pred))
  }
}


mean.cv.errors = apply(cv.errors ,2, mean)
plot(mean.cv.errors,type="b") ; title("Cross Validation PCR")
minM=which(mean.cv.errors==min(mean.cv.errors)) #29 
minM

#--- 3.3 PCR test performance on testdata (with optimal M, on full data)
PcaX <-prcomp(data[,-1])$x 
DataPCR <- data.frame(cbind("meta"=data[,"meta"],PcaX) ) #the full transformed data
f <- as.formula(paste("meta~",paste('PC', 1:minM,sep="",collapse = "+")))
lr.mod.pca <- glm(f,data = DataPCR[train,],family = 'binomial')
PCR.pred= predict(lr.mod.pca, DataPCR[test,],type="response")
PCR.pred<-ifelse(PCR.pred > 0.5 ,1,0)
PCR.misClasificError = mean( (DataPCR[test,"meta"]!= PCR.pred)) 


#--- 3.4 PCR fit on all data with optimal M and get coefficients
PcaXFinal <-prcomp(data[,-1])$x
DataPCR <- data.frame(cbind("meta"=data[,"meta"],PcaXFinal)  )

f <- as.formula(paste("meta~",paste('PC', 1:minM,sep="",collapse = "+")))
lr.mod.pca.final <- glm(f,data = DataPCR,family = 'binomial')
summary(lr.mod.pca.final)
plot(1:29,lr.mod.pca.final$coefficients[-1]);title("coefficients PCR fit")


#--- 3.5 PCR plot performance ROC
pr <- prediction(matrix(PCR.pred), DataPCR[test,1])
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf);title("ROC curve PCR Fit")  #plot ROC curve
auc.pcr <- performance(pr, measure = "auc")
auc.pcr <- auc.pcr@y.values[[1]]  # area under curve
auc.pcr

#summary
lasso.misClasificError
length(lasso.coef[lasso.coef!=0]) 
auc.lasso
ridge.misClasificError
auc.ridge
PCR.misClasificError
minM
auc.pcr


# 4. Examine the effect of different folds on ridge and lasso 
#--------------------------------------------------------------------------
#--- 4.1 make a function for lasso and ridge with folds as a variable
bestlam <- lasso.mod$lambda.min
lasso.coef= predict (lasso.mod.full, s=bestlam, type ="coefficients")
lassoGenes10Folds <- names(data[which(lasso.coef!=0)])[-1]  #get genes of default 10 folds for comparison


glmnet_iter_folds<- function (nfolds, alpha){
  set.seed (1) 
  glmnet.mod =cv.glmnet(x[train,], y[train], alpha =alpha, nfolds=nfolds, family="binomial")
  bestlam <- glmnet.mod$lambda.min
  glmnet.mod.full <- glmnet(x,y, alpha =alpha,family="binomial") 
  glmnet.coef= predict(glmnet.mod.full, s=bestlam, type ="coefficients")
  
  if (alpha==1){ # for lasso we take every non zero gene
  genes <- names(data[which(glmnet.coef!=0)])[-1]
  }
  else { # for ridge we only take the 20 genes with the highest coefficients
    ridge.coef.sort <- data.frame(genes = names(data[,-1]),coefficients=abs(glmnet.coef[-1]));
    genes<- as.vector(ridge.coef.sort[order(-ridge.coef.sort$coefficients),][1:20,1])
  }
  nbSimGenes <-0
  nbDiffGenes <-0
  SimGenes= rep(FALSE,length(genes))
  for (i in seq(1,length(genes))){ 
    flagFoundGene = FALSE;
    for (j in seq(1,length(lassoGenes10Folds))){
      if (genes[i]==lassoGenes10Folds[j]) 
        {nbSimGenes<-nbSimGenes + 1;
        flagFoundGene = TRUE;
        SimGenes[i]=TRUE;
        break} 
    }
    if(!flagFoundGene) {nbDiffGenes<-nbDiffGenes + 1} 
  }
  
  list(nbSimGenes=nbSimGenes,nbDiffGenes=nbDiffGenes,genes=genes,SimGenes = genes[SimGenes])
}



#--- 4.2 iterate over a number of folds

nbFoldsTest <-20 # I will test 20 different folds
nbFoldsComLasso= matrix(ncol=3,data=rep(0,nbFoldsTest*3)) #initialise
colnames(nbFoldsComLasso) = c('folds','nbSimGenes','nbDiffGenes')
rownames(nbFoldsComLasso) = paste(3:(nbFoldsTest+2),c("folds"),  sep = " ")
for (folds in seq(3,2+nbFoldsTest)){
  glmnetFit<-glmnet_iter_folds(folds,1)
  nbFoldsComLasso[folds-2,1] <- folds
  nbFoldsComLasso[folds-2,2] <- glmnetFit$nbSimGenes
  nbFoldsComLasso[folds-2,3] <- glmnetFit$nbDiffGenes
}

glmnetFit<-glmnet_iter_folds(10,0)
glmnetFit$SimGenes

#--- 4.3 look at top 20 genes of Ridge and compare with lasso
nbFoldsTest <-20 # different folds
nbFoldsComRidge= matrix(ncol=3,data=rep(0,nbFoldsTest*3)) #initialise
colnames(nbFoldsComRidge) = c('folds','nbSimGenes','nbDiffGenes')
rownames(nbFoldsComRidge) = paste(3:(nbFoldsTest+2),c("folds"),  sep = " ")
for (folds in seq(3,2+nbFoldsTest)){
  glmnetFit<-glmnet_iter_folds(folds,0)
  nbFoldsComRidge[folds-2,1] <- folds
  nbFoldsComRidge[folds-2,2] <- glmnetFit$nbSimGenes
  nbFoldsComRidge[folds-2,3] <- glmnetFit$nbDiffGenes
}
# lines(x=nbFoldsComRidge[,"folds"],y=nbFoldsComRidge[,"nbDiffGenes"],col="red",lty=2)


#--- 4.4 compare lasso with the default (=10folds) and plot
plot(x=nbFoldsComLasso[,"folds"],y=nbFoldsComLasso[,"nbSimGenes"],type="l",col="green",ylim=c(0,15),lwd = 2,
     xlab="number of folds used in CV",ylab="number of genes")
lines(x=nbFoldsComLasso[,"folds"],y=nbFoldsComLasso[,"nbDiffGenes"],col="red")
lines(x=nbFoldsComRidge[,"folds"],y=nbFoldsComRidge[,"nbSimGenes"],col="green",lty=2)
legend("topright",c('lasso similar', 'lasso different', 'ridge similar(top20'),  lty=c(1,1,2),col=c("green", "red","green"),cex=0.7)
title(main='Similarity in number of genes with lasso 10 fold CV')



# 5. Examine the high correlation genes
#--------------------------------------------------------------------------
HighCorPairs<-which(corMat > 0.9, arr.ind = T)
HighCorPairs<-HighCorPairs[HighCorPairs[,1]<HighCorPairs[,2],]  #filter out diagonals and mirrorvalues
HighCorPairs<-HighCorPairs + 1 #match indexes of Data and coefficient vector, so add 1

HighCorDf <- data.frame(idx1=HighCorPairs[,1], idx2=HighCorPairs[,2])
HighCorDf$Lasso1 <- lasso.coef[HighCorDf$idx1]
HighCorDf$Lasso2 <- lasso.coef[HighCorDf$idx2]
HighCorDf$Ridge1 <- ridge.coef[HighCorDf$idx1]
HighCorDf$Ridge2 <- ridge.coef[HighCorDf$idx2]

sum(HighCorDf$Lasso1>0)
sum(HighCorDf$Lasso2>0)
plot(HighCorDf$Ridge1,HighCorDf$Ridge2,xlab="coeff gene ",ylab=("coeff correlated gene"))
title("Ridge regression : compare coefficients of highly correlated genes",cex=0.5)



