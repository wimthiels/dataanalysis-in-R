
install.packages("ggplot2")
install.packages("lme4")
install.packages("arm")
install.packages("pbkrtest")

library(ggplot2)
library(lme4)
library(lattice)
library(arm)
library(car)
library(pbkrtest)
library(stats)

sleep<-read.delim("sleep.txt", header = TRUE, sep = "\t",
                  dec = ".", fill = TRUE)
attach(sleep)

# EXPLORATORY ANALYSIS
###########################################################################


## Spaghettiplot
  par(mfrow=c(1,1))
  interaction.plot(Days,Subject,Reaction, xlab="Days", ylab="Reaction time", legend=F)
  mean_reactiontimes<-aggregate(Reaction, by=list(sleep$Days), FUN=mean)
  mean_reactiontimes$Group.1 <-mean_reactiontimes$Group.1+1
  lines(mean_reactiontimes,col='red',lwd=2)


# OLS for every trucker
  #sleepplot<- ggplot(sleep, aes(x=Days, y=Reaction, colour=factor(Subject))) + geom_point()
  sleepplot<- ggplot(subset(sleep,Subject<999), 
                   aes(x=Days, y=Reaction, factor(Subject))) +
            geom_point()
  sleepplot2<-sleepplot +  facet_wrap(~Subject,ncol=6) + stat_smooth(method="lm")
  sleepplot2

## creating histograms for the OLS per trucker
  ## Coefficients
  lin.reg.coef <- by(sleep, sleep$Subject,
                     function(data) coef(lm(Reaction ~ Days, data=data)))
  lin.reg.coef1 <- unlist(lin.reg.coef)
  names(lin.reg.coef1) <- NULL
  lin.reg.coef2=matrix(lin.reg.coef1,length(lin.reg.coef1)/2,2,byrow = TRUE)


  ## R squared
  lin.reg.r.squared <- by(sleep, sleep$Subject,
                        function(data) summary(lm(Reaction ~ Days, data=data))$r.squared )
  lin.reg.r.squared1<- as.vector(unlist(lin.reg.r.squared))

   ## Histograms
   par(mfrow=c(2,1))
  hist(lin.reg.coef2[,1],xlab="Intercept",col="lightblue",main="Histogram of individual intercepts")
  hist(lin.reg.coef2[,2],xlab="Slope",col="lightblue",main="Histogram of individual slopes")
  par(mfrow=c(1,1))
  hist(lin.reg.r.squared1,xlab="R squared",col="lightblue",main="Histogram of individual R squared")
  par(mfrow=c(1,1))

  
  ## standardized residuals for within subject OLS
  lin.reg.stan.res <- by(sleep, sleep$Subject,
                     function(data) {lmX<- summary(lm(Reaction ~ Days, data=data) )
                       (lmX$residuals/lmX$sigma)})
  lin.reg.stan.res1 <- unlist(lin.reg.stan.res)
  names(lin.reg.stan.res1) <- NULL
  lin.reg.stan.res2=matrix(lin.reg.stan.res1,length(lin.reg.stan.res1)/10,10,byrow = TRUE)

  plot(0:9,lin.reg.stan.res2[1,],pch=" ",
       main='within subject OLS residuals (standardized)',
       ylim=c(-2.5,2.5),
       xlab='Days',ylab='residual (standardized)',
       cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=1.5)
  
  
  for (i in 1:18) {
    lines(0:9,lin.reg.stan.res2[i,])
  }
  stan_error_mean= rep(0,10)
  for (i in 1:10) {
    stan_error_mean[i]=mean(lin.reg.stan.res2[,i])
  }
  lines(0:9,stan_error_mean,col='red',lwd=2)


#making a variance plot (based on residuals from the individual regression lines)
  trucker_res_error <- data.frame(matrix(data = 0,nrow = 1,ncol = 10))
  colnames(trucker_res_error) <- c("day0","day1","day2","day3","day4","day5","day6","day7","day8","day9")
  #collect all the residuals
  trucker_counter <- 0
  for (i in unique(Subject)) {
    trucker_counter <- trucker_counter + 1 
    lmTrucker<-lm(Reaction~Days,data=subset(sleep,Subject==i))
    trucker_res_error[trucker_counter,]<-lmTrucker$residuals
    rownames(trucker_res_error)[trucker_counter]<-c(i)
  }
   #calculate variance for every day
  res_variance_per_day <- c(rep(0,10))
  for (i in 1:10) {
    res_variance_per_day[i] <- var(trucker_res_error[,i])
  }
  trucker_res_error[,7]
  var(trucker_res_error[,7])
  plot(0:9,res_variance_per_day,
     xlab='day number',ylab='residual variance (after OLS per trucker)',
     cex=0.8,cex.lab=0.9, cex.axis=0.8, cex.main=1.5)
  lines(0:9,res_variance_per_day)

  
  ## histogram for the st dev (sigma) of each OLS : large variation ?
  lin.reg.sigma <- by(sleep, sleep$Subject,
                         function(data) {lmX<- summary(lm(Reaction ~ Days, data=data) )
                         (lmX$sigma)})
  par(mfrow=c(1,1))
  hist(as.matrix(lin.reg.sigma),xlab="residual error standard deviation",col="lightblue",
       main="Individual residual error standard deviation")


  

# FITTING THE DATA
###########################################################################

  sleep.lmer<-lmer(Reaction~1+Days +(1 + Days|Subject), REML = TRUE,data=sleep)
  summary(sleep.lmer)
  
  confint(sleep.lmer,par=5:6,method="Wald",oldNames = FALSE)
  anova(sleep.lmer)
  
  

# INFERENCE
###########################################################################
## inference for the fixed effects  //////
  #Get the KR-approximated degrees of freedom
  sleep.lmer.df.KR <- get_ddf_Lb(sleep.lmer, fixef(sleep.lmer))

  #Get p-values from the t-distribution using the t-values and approximated degrees of freedom
  sleep.lmer.coef=coef(summary(sleep.lmer))
  sleep.lmer.p.KR <- cbind(sleep.lmer.coef,'p-value'=2 * (1 - pt(abs(sleep.lmer.coef[,3]), sleep.lmer.df.KR)))
  sleep.lmer.p.KR

## inference for the variance effects////////////
  # leave out the covariance ?

  sleep.lmer<-lmer(Reaction~1+Days +(1 + Days|Subject), REML = TRUE, data=sleep)
  sleep.lmer.null <- lmer(Reaction ~ 1 + Days + (1|Subject) + (0+Days|Subject), data=sleep, REML=TRUE)
  anova(sleep.lmer,sleep.lmer.null,test=TRUE)

  # leave out the variance of the slope ?
 
  sleep.lmer      <- lmer(Reaction~1+Days +(1 + Days|Subject), REML = TRUE, data=sleep)
  sleep.lmer.null <- lmer(Reaction ~ 1 + Days + (1|Subject), REML=TRUE, data=sleep)
  logratio_statistic <- -2 * (logLik(sleep.lmer.null,REML=TRUE) - logLik(sleep.lmer,REML=TRUE))
  P_chi_df1 <-as.numeric(pchisq(logratio_statistic, df=1,lower.tail = FALSE))
  P_chi_df2 <-as.numeric(pchisq(logratio_statistic, df=2,lower.tail = FALSE))
  P_null_model <- 1/2 * (P_chi_df1 + P_chi_df2)
  P_null_model

  # the FINAL MODEL
  ###########################################################################
## State the FINAL MODEL (so after inferencetests)
  sleep.lmer.final <- lmer(Reaction ~ 1 + Days + (1|Subject) + (0+Days|Subject), data=sleep, REML=TRUE)
  summary(sleep.lmer.final)
  
  
  sleep.lmer.final.VarCorr<-VarCorr(sleep.lmer.final)
  sleep.lmer.final.VarCorr[]
  
  confint(sleep.lmer.final,par=3:5,method="Wald",oldNames = FALSE) ## equals 251.405-(1.96*6.885)
  
  # INFERENCE (part deux => random effects)
  ###########################################################################
  
# Inference for the random effects/////

  ## Predicted random effects (only the b's , so without fixed effect)
  
  sleep.lmer.re=ranef(sleep.lmer.final)$Subject
  sleep.lmer.re 

  # show the shrinkage effect
  # -------------------------
    # LMM : get the b's from our multilevel model
  ind.coef=coef(sleep.lmer.final)$Subject # the b's + fixed effects
  int.subject=ind.coef[,1]
  slope.subject=ind.coef[,2]
  

     #OLS : get the regression coefficients of the individual linear regressions
  trucker_reg_coeff <- data.frame(rep(0,18),0)
  colnames(trucker_reg_coeff) <- c("intercept","Days")
  trucker_counter <- 0
  for (i in unique(Subject)) {
    trucker_counter <- trucker_counter + 1 
    lmTrucker<-lm(Reaction~Days,data=subset(sleep,Subject==i))
    trucker_reg_coeff[trucker_counter,]<-lmTrucker$coefficients
    rownames(trucker_reg_coeff)[trucker_counter]<-c(i)
  }

    #plotting the shrinkage effect
  plot(trucker_reg_coeff[,2],trucker_reg_coeff[,1],col='red',
     main="per subject OLS versus LMM estimates of the slope and intercepts",
     cex.lab=0.8, cex.axis=0.9, cex.main=0.7, cex.sub=0.9,
     pch=1,
     xlab='slope',ylab='intercept')
  points(ind.coef[,2],ind.coef[,1],pch=19)

  for (i in seq(1,18)) {
  arrows(trucker_reg_coeff[i,2],trucker_reg_coeff[i,1],
         ind.coef[i,2],ind.coef[i,1], col= 'red',length=0.1)
  }                    

  legend("topright", c('OLS', 'LMM','population'), pch=c(1,19,13),cex=0.5,
     col=c('red','black','blue'),xpd=TRUE)

  lm_plain<-lm(Reaction~Days,data=sleep)$coefficients
  points(lm_plain[2],lm_plain[1],col='blue',pch=13)

  truckersID<-rownames(trucker_reg_coeff)
  truckersID<-ifelse(truckersID=='335',"335","")
  text(trucker_reg_coeff[,2], trucker_reg_coeff[,1], labels=truckersID,cex=0.5,pos=3)

  ## sanity check : normality of the residual errors after fitting the model
  # get the estimated mean values from our model : y-hat

  
  #get standardised final residuals
  final_resid<-summary(sleep.lmer.final)$residuals
  
  # just a check to see that these are indeed the standardized res = ok 
  cov_matrix <- rbind (rep(1,10),0:9)
  Y_hat_model <- as.matrix(ind.coef)  %*%  cov_matrix 
  
  par(mfrow=c(1,1))
  plot(0:9,Y_hat_model[1,],pch="",ylim=cbind(0,500))
  for (i in 1:18) lines(0:9,Y_hat_model[i,])
  
  sleep_matrix<-matrix(sleep$Reaction,ncol=10,byrow=TRUE)
  final_resid2 <- (-Y_hat_model + sleep_matrix)/summary(sleep.lmer.final)$sigma
  final_resid2 ## equals final_resid : so ok !
  
  
  final_resid2
  outlier_points<- ifelse(abs(final_resid2)>3,"outlier","ok")
  outlier_points
  qqnorm(final_resid2, main='QQ plot residuals final model',ylab='standardised residuals')
  qqline(final_resid2, col = 2) #The same 3 outlier datapoints: ok


  
  
  
  
  
  
  
  
