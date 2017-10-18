#0# preparation#####

#install.packages("GGally")
#install.packages("gam")
#install.packages("ROCR")
#install.packages('pROC')
library(GGally)
library(ggplot2)
library(boot)
library(gam)
library(ROCR)
library(pROC)

setwd(dir="C:/Users/Administrator/Documents/OneDrive/R_MyFiles")
heart<-read.csv("SAheart.data")
set.seed(1)
par ( mfrow =c(1,1) )

#1# explore the data####


heart<-heart[,-1]  #get rid of row.names column
attach(heart)
dim(heart)
names(heart) #chd is myocardial infarction = response = binary
sum(is.na(heart))   # no missing data
summary(heart)
plot(heart)
ggpairs(heart)
sum(heart$chd==0)


#2# looking at family history as predictor ####

lrmod1<-glm(chd~famhist,data=heart,family = 'binomial')
coefFamhist<-lrmod1$coefficients[[2]]
Intercept<-lrmod1$coefficients[[1]]
RiseInP<-(exp(coefFamhist+Intercept)/(1+exp(coefFamhist+Intercept))) - 
  (exp(Intercept)/(1+exp(Intercept)))
RiseInP


#3# glm with backward selection      ####
#--step0 : functions-------------------------------------------------
get_formula <-function(predselect,modeltype){
  
  if (modeltype=="glm"){
    predictors = c("1",colnames(heart)[predselect])
    formula<-paste("chd~",paste(predictors, collapse = "+"),sep = "")
  }
  
  
  if (modeltype=="gam"){
    
    famHistIdx = which(names(predsel_df)=="famhist")
    famHistSwitch = predselect[famHistIdx]
    famHistSwitch
    famHistIdx
    predselect[famHistIdx] = FALSE  #turn of famhist, not via spline
    predictorsSpline = c(colnames(heart)[predselect])
    predselect[famHistIdx] = TRUE
    
    formula  <-"chd~1"
    
    if (famHistSwitch){
      formula <- paste(formula,"+famhist")
    }
    
    if (length(predictorsSpline) != 0) {
      formula<-  paste(formula,"+",
                       paste("s(",predictorsSpline, ",4)"
                             ,collapse = "+",sep="")
      )
      
    }
  }
  as.formula(formula)
}

determine_pred_2Bdropped <- function(predselect,modeltype){  #find predictor that can be removed from the model with the least increase in the deviance as a result
  
  minDevianceIdx = 0
  minDeviance = 999999999999999
  
  for (j in seq(1,length(predselect))) {
    if (predselect[j]==TRUE){
      
      predselect[j]=FALSE  # switch off this predictor and see how good the fit is using CV
      
      if(modeltype=="glm"){
        modFit <-glm(get_formula(predselect,"glm"),data=heart)}
      if(modeltype=="gam"){
        modFit <- gam(get_formula(predselect,"gam"),data=heart)}
      
      ###
      
      
      if (modFit$deviance < minDeviance) { #keep track of the best model
        minDevianceIdx=j
        minDeviance = modFit$deviance 
      }
      
      predselect[j]=TRUE  # switch'm back on
    }
  }
  
  list(predDropIdx = minDevianceIdx, remainingDeviance = minDeviance )
}


#--step3.1 glm: set up the datastructures-----------------------------------------
predselect <- rep(x=TRUE,9)
predselect[10] <- FALSE
names(predselect)<-colnames(heart) 
predselect  #a vector to indicate columns that are in the model

cv.10lrmods <- rep(x=0,10)
names(cv.10lrmods)<-paste(seq(9,0,-1),"pred",sep="") 
cv.10lrmods  #a vector to store the CV errors for every model

predSelHisGlm <- data.frame(matrix(data = FALSE,nrow = 10,ncol = 10),row.names = names(cv.10lrmods))
colnames(predSelHisGlm) <-colnames(heart)
predSelHisGlm   #a matrix to keep track of the predictors selected in every model


#--step3.2 glm: do backward selection, every iteration loosing a predictor-----------------------------
set.seed(1)
for (i in seq(1,10)){  
  
  predSelHisGlm[i,]<-predselect
  
  # fit the logistic regression and get the cv-error to able to select optimal nb of pred.
  lrmod    <-glm(get_formula(predselect,"glm"),data=heart)
  
  cv.10lrmods[i]<-cv.glm(data=heart, 
                         glmfit=lrmod,
                         cost=function(trueY, predY = 0) mean(abs(trueY-predY) > 0.5),
                         K=10)$delta[1]
  
  
  if (i!=10){
    predDropIdx = determine_pred_2Bdropped(predselect,"glm") 
    predselect[predDropIdx[[1]]] = FALSE; #loose that predictor
  }
  
}

#--step3.3 glm: plot the CV of the backward selection and get best nb predictors-------------------------

plotpredSelHisGlm<- ifelse(predSelHisGlm==TRUE,"*"," ")
plotpredSelHisGlm[-10,-10] # show overview matrix

plot(x=seq(9,0,-1),y=cv.10lrmods,
     main="CV error GLM backwards selection",
     xlab="nb of predictors in model",
     ylab="CV error (=misclassification %)",
     xlim=c(9,0)
)
which.min(cv.10lrmods)
points(5, cv.10lrmods[5], col = "red",pch=3,lwd=5)
lines(x=seq(9,0,-1),y=cv.10lrmods)

get_formula(as.logical(predSelHisGlm["5pred",]),"glm")

glmFitFull <- glm(get_formula(as.logical(predSelHisGlm["5pred",]),"glm"),data=heart)
summary(glmFitFull)

#--step3.4 glm : explain the variables that were dropped from the model : why ?----------------------
#------a : examine correlation with response----
predsel<-as.logical(predSelHisGlm["5pred",])
predsel_df<-predSelHisGlm["5pred",]
predsel_df
corWithRes <- matrix(rep(0,10),nrow=1,dimnames = list("correlation",colnames(heart)))
corWithRes[,predsel][-3] <- cor(heart[,"chd"],heart[,predsel][,-3]) #exclude famhist
corWithRes[,!predsel] <- cor(heart[,"chd"],heart[,!predsel])
corWithRes[,"famhist"] <- NA
corWithResDf<-as.data.frame(t(corWithRes))
corWithResDf

plotData <-subset(corWithResDf,subset=!rownames(corWithResDf) %in% c('chd','famhist'))
plotPred <-as.logical(subset(predsel_df, select=-c(chd,famhist)))
plot<-ggplot(plotData,aes("",correlation,colour=plotPred)) 
plot<- plot + geom_jitter(size=1,width=0.03)
plot <- plot + labs(title = "Compare correlation of predictors with response(chd)", 
                    x = "", y = "Correlation with chd", color = "Pred in Model ?")
plot <- plot + geom_text(aes(label = rownames(plotData)),check_overlap = FALSE,nudge_x = -0.1)
plot <- plot + scale_color_manual(labels = c("NO", "YES"), values = c("red", "blue"))
plot <- plot +theme_bw()

plot

#------b : examine correlation with adiposity----

predsel<-as.logical(predSelHisGlm["8pred",])
predsel_df<-predSelHisGlm["8pred",]
predsel_df
corWithAdip <- matrix(rep(0,10),nrow=1,dimnames = list("correlation",colnames(heart)))
corWithAdip[,predsel][-5] <- cor(heart[,"adiposity"],heart[,predsel][,-5]) #exclude famhist
corWithAdip[,!predsel] <- cor(heart[,"adiposity"],heart[,!predsel])
corWithAdip[,"famhist"] <- NA
corWithAdipDf<-as.data.frame(t(corWithAdip))
corWithAdipDf

predsel<-as.logical(predSelHisGlm["5pred",])
predsel_df<-predSelHisGlm["5pred",]
plotData <-subset(corWithAdipDf,subset=!rownames(corWithAdipDf) %in% c('adiposity','famhist'))
plotPred <-as.logical(subset(predsel_df, select=-c(adiposity,famhist)))
plot<-ggplot(plotData,aes("",correlation,colour=plotPred)) 
plot<- plot + geom_jitter(size=1,width=0.03)
plot <- plot + labs(title = "Compare correlation of predictors with adiposity", 
                    x = "", y = "Correlation with adiposity", color = "Pred in Model ?")
plot <- plot + geom_text(aes(label = rownames(plotData)),check_overlap = FALSE,nudge_x = -0.1)
plot <- plot + scale_color_manual(labels = c("NO", "YES"), values = c("red", "blue"))
plot <- plot +theme_bw()

plot

#------c : examine correlation with sbp----

predsel<-as.logical(predSelHisGlm["7pred",])
predsel_df<-predSelHisGlm["7pred",]

corWithSbp <- matrix(rep(0,10),nrow=1,dimnames = list("correlation",colnames(heart)))
corWithSbp[,predsel][-5] <- cor(heart[,"sbp"],heart[,predsel][,-4]) #exclude famhist
corWithSbp[,!predsel] <- cor(heart[,"sbp"],heart[,!predsel])
corWithSbp[,"famhist"] <- NA
corWithSbpDf<-as.data.frame(t(corWithSbp))
corWithSbpDf

predsel<-as.logical(predSelHisGlm["5pred",])
predsel_df<-predSelHisGlm["5pred",]
plotData <-subset(corWithSbpDf,subset=!rownames(corWithSbpDf) %in% c('sbp','famhist'))
plotPred <-as.logical(subset(predsel_df, select=-c(sbp,famhist)))
plot<-ggplot(plotData,aes("",correlation,colour=plotPred)) 
plot<- plot + geom_jitter(size=1,width=0.03)
plot <- plot + labs(title = "Compare correlation of predictors with systolic blood pressure", 
                    x = "", y = "Correlation with systolic blood pressure", color = "Pred in Model ?")
plot <- plot + geom_text(aes(label = rownames(plotData)),check_overlap = FALSE,nudge_x = -0.1)
plot <- plot + scale_color_manual(labels = c("NO", "YES"), values = c("red", "blue"))
plot <- plot +theme_bw()

plot

#4# GAM with backward selection      ####
#--step4.1 gam: initialise the datastructures #-----------
predselect <- rep(x=TRUE,9)
predselect[10] <- FALSE

aic.10gams <- rep(x=0,10)
names(aic.10gams)<-paste(seq(9,0,-1),"pred",sep="") 
aic.10gams  #a vector to store the AIC for every model

predSelHisGam <- data.frame(matrix(data = FALSE,nrow = 10,ncol = 10),row.names = names(cv.10lrmods))
colnames(predSelHisGam) <-colnames(heart)


#--step4.2 gam : do backward selection, every iteration loosing a predictor----------
set.seed(1)
for (i in seq(1,10)){  
  
  print(paste("i->",i))
  predSelHisGam[i,]<-predselect
  
  # fit the logistic regression and get the cv-error to able to select optimal nb of pred.
  gamMod    <-gam(get_formula(predselect,"gam"),data=heart)
  aic.10gams[i]<-gamMod$aic
  
  if (i!=10){
    predDropIdx = determine_pred_2Bdropped(predselect,"gam") 
    predselect[predDropIdx[[1]]] = FALSE; #loose that predictor
  }
  
}

#--step4.3 gam: plot the AIC of the backward selection and get best nb predictors-------------------------

plotpredSelHisGam<- ifelse(predSelHisGam==TRUE,"*"," ")
plotpredSelHisGam[-10,-10] # show overview matrix
plot(predSelHisGam)

plot(x=seq(9,0,-1),y=aic.10gams,
     main="AIC GAM backwards selection",
     xlab="nb of predictors in model",
     ylab="AIC",
     xlim=c(9,0)
)
which.min(aic.10gams)
points(6, aic.10gams[which.min(aic.10gams)], col = "red",pch=3,lwd=5)
lines(x=seq(9,0,-1),y=aic.10gams)

get_formula(as.logical(predSelHisGam["6pred",]),"gam")

#--step4.4 gam: check for non-linearity--------------
gamFitFull <- gam(get_formula(as.logical(predSelHisGam["6pred",]),"gam"),data=heart)
summary(gamFitFull)

get_formula(as.logical(predSelHisGam["5pred",]),"gam")
#loose obesity, to get the final gam model
gamFitFull <- gam(get_formula(as.logical(predSelHisGam["5pred",]),"gam"),data=heart)

#5# Compare performance of gam and glm     ####
#--- step 5.1 fit gam and glm on training data ----
train = sample(x=nrow(heart),size = (2*nrow(heart))/3,replace=FALSE)


glmFit <- glm(get_formula(as.logical(predSelHisGlm["5pred",]),"glm"),data=heart[train,])
gamFit <- gam(get_formula(as.logical(predSelHisGam["6pred",]),"gam"),data=heart[train,])
summary(glmFit)
summary(gamFit)

#--- step 5.2 compare  the misclassification rate ----
table(heart[-train,][,"chd"])
dim(heart[-train,])
glmPredict <- predict(glmFit,newdata=heart[-train,],type="response")
glmPredict <- ifelse(glmPredict>0.5,1,0)
misClasificErrorGlm <- mean(glmPredict != heart[-train,][,"chd"])
print(paste('Accuracy GLM:',1-misClasificErrorGlm))
print(paste('Misclassification rate GLM :',misClasificErrorGlm))

gamPredict <- predict(gamFit,newdata=heart[-train,],type="response")
gamPredict <- ifelse(gamPredict>0.5,1,0)
misClasificErrorGam <- mean(gamPredict != heart[-train,][,"chd"])
print(paste('Accuracy GAM :',1-misClasificErrorGam))
print(paste('Misclassification rate GAM :',misClasificErrorGam))


#--- step 5.3 plot ROC curves ----
pr<-prediction(as.numeric(glmPredict),heart[-train,][,"chd"])
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
title("ROC curve glm")  #plot ROC curve
auc.glm <- performance(pr, measure = "auc")
auc.glm <- auc.glm@y.values[[1]]  # area under curve
print(paste('AUC GLM :',auc.glm))
#other package
rocGlm=roc(heart[-train,][,"chd"], as.numeric(glmPredict))
plot(rocGlm)


pr<-prediction(as.numeric(gamPredict),heart[-train,][,"chd"])
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
title("ROC curve gam")  #plot ROC curve
auc.gam <- performance(pr, measure = "auc")
auc.gam <- auc.gam@y.values[[1]]  # area under curve
print(paste('AUC GAM :',auc.gam))


#--- step 5.4 the final model  ----
gamFinal<- gam(chd ~ 1 + tobacco + ldl + famhist + s(typea,4) + age,data=heart)

gamFinalPredict <- predict(gamFinal,newdata=heart[-train,],type="response")
gamFinalPredict <- ifelse(gamFinalPredict>0.5,1,0)
misClasificErrorGamFinal <- mean(gamFinalPredict != heart[-train,][,"chd"])
print(paste('Accuracy final model:',1-misClasificErrorGamFinal))
print(paste('Misclassification rate final model :',misClasificErrorGamFinal))

pr<-prediction(as.numeric(gamFinalPredict),heart[-train,][,"chd"])
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
title("ROC curve gamFinal")  #plot ROC curve
auc.gamFinal <- performance(pr, measure = "auc")
auc.gamFinal <- auc.gamFinal@y.values[[1]]  # area under curve
print(paste('AUC gamFinal :',auc.gamFinal))


#--- step 5.5 compare the backfitting procedure in gam and glm ----

plotBarBackFit <- data.frame(modeltype=character(),
                             predictor=character(),stringsAsFactors=FALSE)

for (row in row.names(predSelHisGlm)){
  for (column in names(predSelHisGlm)){
    if(predSelHisGlm[row,column]){
      newrow = c("GLM",column)
      plotBarBackFit[nrow(plotBarBackFit)+1,] <-newrow
    }}
}
for (row in row.names(predSelHisGam)){
  for (column in names(predSelHisGam)){
    if(predSelHisGam[row,column]){
      newrow = c("GAM",column)
      plotBarBackFit[nrow(plotBarBackFit)+1,] <-newrow
    }}
} 


barChart<-ggplot(plotBarBackFit, aes(predictor, fill = modeltype))
barChart<-barChart + geom_bar(position = "dodge")
barChart<- barChart +theme_bw()
barChart <- barChart + labs(title = "compare the inclusion of predictors during backward selection", 
                            x = "Predictor", y = "included in how many models ?")

barChart<-barChart + theme(
  # axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())
barChart


#--- step 5.6 do an anova test on all models ----
anova (glmFitFull ,gamFinal,gamFitFull, test ="F")
glmFitFull$aic
gamFinal$aic
gamFitFull$aic

#--- step 5.7 plotting all models ----
par ( mfrow =c(1,1) )
par ( mfrow =c(2,3) )
plot.gam(glmFitFull,se=T,mfrow=c(2,3),col="red")
par ( mfrow =c(1,1) )
par ( mfrow =c(2,3) )
plot.gam(gamFitFull,se=T,mfrow=c(2,3),col="red")
par ( mfrow =c(1,1) )
par ( mfrow =c(2,3) )
plot.gam(gamFinal,se=T,mfrow=c(2,3),col="red")



