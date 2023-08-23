## ---------------------------------------------------------------------------------------------------------
################################################################################
# Required packages
################################################################################
rm(list=ls())
library(mombf)
library(survival)
library(xtable)
library(ggplot2)
library(BhGLM)
library(BVSNLP)

################################################################################
# Data preparation
################################################################################
dat = read.table("ftbrs_survdata.txt",header=TRUE)
y= Surv(log(dat$Timerecurrence/12), event=dat$recurrence)
X= dat[,-1:-3]; X$stage= factor(X$stage)
X[,-1] <- apply(X[,-1],2,scale)
df <- data.frame(y = y, as.data.frame(X), time = dat$Timerecurrence/12, status = dat$recurrence)


# Visualise the Kaplan-Meier estimator (of survival)
km <- survfit(Surv(time, status) ~ 1, data = df)
plot(km, xlab = "Years", ylab = "Survival", main = "Kaplan-Meier estimator",
     cex.axis = 1.5, cex.lab = 1.5)


## ---------------------------------------------------------------------------------------------------------
####################################################################################################
##  DATA ANALYSIS FOR MODEL ONLY WITH TGFB AND STAGE
## TGFB is supposed to be an important variable in cancer tumour evolution
####################################################################################################

#Run Bayesian variable selection

ms0= modelSelection(y ~ X$stage + X$tgfb,  
                    priorCoef= momprior(taustd = 1), 
                    priorDelta= modelbbprior(1,1))
pp0= postProb(ms0)
pp0
ms0$margpp #TGFB gets 0.9014 inclusion probability


## ---------------------------------------------------------------------------------------------------------
################################################################################
# Bayesian Variable Selection (BVS) for TGFB and other variables
################################################################################

# Regression formula with linear effects
f1y <-  formula(paste('y ~ ',paste(colnames(X), collapse='+'),sep=''))

#-------------------------------------------------------------------------------
# BVS: Prior 1 - Zellner's prior
#-------------------------------------------------------------------------------
# Priors for the coefficients
priorCoef <- zellnerprior(tau=nrow(X))
# Prior on the model space
priorDelta <- modelbbprior(1,1)
# Selection step
ms1 <- modelSelection(f1y,data=df,
                      priorCoef=priorCoef,priorDelta=priorDelta,
                      enumerate=FALSE, method='Laplace', niter=10000, 
                      center = TRUE, scale = TRUE)

# Calculating model posterior probabilities
pp1 <- postProb(ms1)
# Top models
head(pp1)
# Marginal inclusion probabilities
mp1 <- ms1$margpp

# Visualising the marginal inclusion probabilities > 0.5
ind1 <- which(mp1>0.5)
df1<-data.frame(names=names(mp1[ind1]),mp = as.vector(round(mp1[ind1],digits = 3)))
df1

ggplot(df1,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")


#-------------------------------------------------------------------------------
# BVS: Prior 2 - MOM prior (Laplace)
#-------------------------------------------------------------------------------
# Priors for the coefficients
priorCoefm <- momprior(taustd=1) # unit information p-mom (or tau = 0.192)
# Prior on the model space
priorDelta <- modelbbprior(1,1)

priorGroup= groupzellnerprior(tau=nrow(X))
# Selection step
ms2 <- modelSelection(f1y,data=df,
                      priorCoef=priorCoefm,priorDelta=priorDelta,
                      priorGroup = priorGroup,
                      enumerate=FALSE, method='Laplace', niter=10000, 
                      center = TRUE, scale = TRUE)



# Calculating model posterior probabilities 
pp2 <- postProb(ms2)
# Top models
head(pp2)
# Marginal inclusion probabilities
mp2 <- ms2$margpp


# Visualising the marginal inclusion probabilities
ind2 <- which(mp2>0.5)
df2<-data.frame(names=names(mp2[ind2]),mp = as.vector(round(mp2[ind2],digits = 3)))
df2

ggplot(df2,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")

# What conclusions can we get from this?
# MoM priors indicate that stage explains recurrence time. Does this suggest any additional 
# strategy in model selection? --- Stratification. 


#-------------------------------------------------------------------------------
# BVS: Prior 3 - MOM prior (ALA)
#-------------------------------------------------------------------------------

priorGroup= groupzellnerprior(tau=nrow(X))
# Selection step
ms3 <- modelSelection(f1y,data=df,
                      priorCoef=priorCoefm,priorDelta=priorDelta,
                      priorGroup = priorGroup,
                      enumerate=FALSE, method='Laplace', niter=10000, 
                      center = TRUE, scale = TRUE)



# Calculating model posterior probabilities 
pp3 <- postProb(ms3)
# Top models
head(pp3)
# Marginal inclusion probabilities
mp3 <- ms3$margpp


# Visualising the marginal inclusion probabilities
ind3 <- which(mp3>0.5)
df3<-data.frame(names=names(mp3[ind3]),mp = as.vector(round(mp3[ind3],digits = 3)))
df3

ggplot(df3,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")



## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Cox PiMOM priors (sensitivity analysis)
#####################################################################

# Survival object with times
y2 <- Surv(time = df$time, df$status)

X2 <- X

bvsfit <- bvs(X = data.frame(X2), resp = as.matrix(cbind(df$time,df$status)), 
              family  = "survival", nlptype = "piMOM", mod_prior = "unif",
              niter = 100)

# Models with highest probability
head(bvsfit$max_models)

# Highest probabilities
head(exp(bvsfit$max_prob_vec - max(bvsfit$max_prob_vec))/sum(exp(bvsfit$max_prob_vec - max(bvsfit$max_prob_vec) )))

# Highest posterior model
bvsfit$HPM

# Estimated coefficients
bvsfit$beta_hat

# Number of Visited Models
bvsfit$num_vis_models

# Visualising the marginal inclusion probabilities
mp = as.vector(round(bvsfit$inc_probs,digits = 3))
dfi <- data.frame( names1 = names(mp1)[-1], mp = as.vector(round(bvsfit$inc_probs,digits = 3)))

index = tail(order(dfi$mp),n=10)
dfint <- dfi[index,]


ggplot(dfint,aes(x =  reorder(names1, -mp), y = mp,label = mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")


## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Spike-and-Slab LASSO (sensitivity analysis)
#####################################################################

# Requires a different design matrix as it does not admit factors
X3= dat[,-1:-3]; 
X3 <- apply(X3,2,scale)

# Stage as a binary variable
#stageM <- matrix(0, ncol = 3, nrow =  dim(X3)[1])
#for(i in 1:dim(X3)[1]) stageM[i,df$stage[i]] = 1

# Adding stage indicators
#X3 <- as.matrix(cbind(stageM[,2:3],X3))

# Selection step
spsl <- bmlasso(x = X3, y =  y2, family = "cox")

# Posterior modes (only option in this package, no MCMC)
head(sort(spsl$coef,decreasing = TRUE))

# Posterior modes
plot.bh(coefs = spsl$coef, threshold = spsl$df, gap = 1, 
        main = "Spike-and-Slab LASSO", lty = 2) 

# Posterior Odds-Ratios for the modes
plot.bh(coefs = spsl$coef, threshold = spsl$df, gap = 1, 
        main = "Spike-and-Slab LASSO", lty = 2, OR = TRUE) 



## ---- warning=FALSE---------------------------------------------------------------------------------------
#####################################################################
# Cox-LASSO (sensitivity analysis)
#####################################################################

# Selection step
cv.fit = try(cv.glmnet(x =  X3, y = y2, family="cox", 
                       maxit=1000, nfolds=10, alpha=1), silent=TRUE)
fit = try(glmnet(x = X3, y=y2, family = "cox", maxit=1e5, alpha=1), silent=TRUE)

# active variables (lambda.min)
b.coxlasso = as.double(coef(fit, s=cv.fit$lambda.min))
which(b.coxlasso!=0)
colnames(X)[which(b.coxlasso!=0)]

# active variables (lambda.1se)
b2.coxlasso = as.double(coef(fit, s=cv.fit$lambda.1se))
which(b2.coxlasso!=0)
colnames(X)[which(b2.coxlasso!=0)]

# Fitting a CoxPH model using the survival R package
fc <- formula(paste('y2 ~ ',paste(colnames(X3)[which(b.coxlasso!=0)], collapse='+'),sep=''))
fit2 <- coxph(fc, data = df)
summary(fit2)



