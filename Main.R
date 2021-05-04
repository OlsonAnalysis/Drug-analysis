#First we extract the data from the file
BigPillData<-get(load("Data-Olson.Rdata"))


# we split the data into males and females in order to see if the same model should be used for both
#Since we will do a PCA analysis, we scale the predictors
BPMale<-subset(BigPillData, IsMale == TRUE)
BPFemale<-subset(BigPillData, IsMale == FALSE)
BPMalePred<-BPMale[,3:6]
BPFemalePred<-BPFemale[,3:6]
scaledMalePred<-scale(BPMalePred)
scaledFemalePred<-scale(BPFemalePred)


#************************************************************************************************************************
#PCA analysis for male predictors
covmat <- cov(scaledMalePred)
egnvectors <- eigen(covmat)$vectors # get eigenvectors
scoresmat <- scaledMalePred %*% egnvectors

egnvalues <-  eigen(covmat)$values
pve1 <- egnvalues / sum(egnvalues) #Each PVE is each eigenvalue over the total sum
cumPVE1 <- c(pve1[1])
for (j in 2:4) {
  cumPVE1 <- c(cumPVE1, cumPVE1[j-1]+pve1[j]) #Incrimental increase in explanation
}
cumPVE1


#PCA analysis for female predictors
covmat <- cov(scaledFemalePred)
egnvectors <- eigen(covmat)$vectors # get eigenvectors
scoresmat <- scaledFemalePred %*% egnvectors

egnvalues <-  eigen(covmat)$values
pve1 <- egnvalues / sum(egnvalues) #Each PVE is each eigenvalue over the total sum
cumPVE1 <- c(pve1[1])
for (j in 2:4) {
  cumPVE1 <- c(cumPVE1, cumPVE1[j-1]+pve1[j]) #Incrimental increase in explanation
}
cumPVE1

#Since the variance explained by each of the variables is roughly the same for both male and female subsets
#this might indicate that the same model would work for both males and females. However, it is possible that although
#the variance explained is similar, some predictors will have a larger effect on the response, a simple linear regression
#should verify if this is the case
#Also, all the predictors explain a large part of the variance so we cannot eliminate any of them.



#**********************************************************************************************************************************************************************************
#We change the True/False values to 0-1 in order to perform various statistical methods
numericBigPill<-BigPillData
cols <- sapply(BigPillData, is.logical)
numericBigPill[,cols] <- lapply(BigPillData[,cols], as.numeric)
head(numericBigPill)


#first let's simply do a linear regression to get a good idea of which factors have a large effect on complications
#let's scale the variables so that the magnitude of the coefficients have meaning

scaleNumBP<-scale(numericBigPill)
lm(scaleNumBP[,1]~scaleNumBP[,2:6])

#all the coefficients magnitude are substantially different than 0 which confirms our PCA analysis
#IE: all the predictors have some effect on the response

#let's see if the magnitude of coefficients are similar for both male and female:
BPnumMale<-subset(numericBigPill, IsMale == 1)
BPnumFemale<-subset(numericBigPill, IsMale == 0)
scaleBPnumMale<-scale(BPnumMale)
scaleBPnumFemale<-scale(BPnumFemale)

lm(scaleBPnumMale[,1]~scaleBPnumMale[,3:6])
lm(scaleBPnumFemale[,1]~scaleBPnumFemale[,3:6])

#the average blood pressure seems to affect the chance of Females getting complications about twice as much as for males
#the neutrophilePercentage seems to have about 10x more effect on females than on males
#the other predictors seem pretty much the same
#this indicates that it might be better to have seperate models for males and females

#to determine if there are interactions between the different variables we compute the correlation

cor(numericBigPill)
cor(BPnumMale[,-2])
cor(BPnumFemale[,-2])

#we note that the correlation between the different predictors are all very small 
#To determine if these are statistically different than 0 we use the psychometric package
#install.packages("psychometric") 
library(psychometric) 

ConfidenceIntervals<-numeric(0)
for(i in 1:5){
  for(j in 1:5){
    ConfidenceIntervals<-c(ConfidenceIntervals, CIr(cor(numericBigBill[,-2])[i,j],n=5000,level=0.95))
  }
}

ConfidenceIntervals

ConfidenceIntervalsMale<-numeric(0)
for(i in 1:5){
  for(j in 1:5){
ConfidenceIntervalsMale<-c(ConfidenceIntervalsMale, CIr(cor(BPnumMale[,-2])[i,j],n=5000,level=0.95))
  }
}

ConfidenceIntervalsMale

ConfidenceIntervalsFemale<-numeric(0)
for(i in 1:5){
  for(j in 1:5){
    ConfidenceIntervalsFemale<-c(ConfidenceIntervalsFemale, CIr(cor(BPnumFemale[,-2])[i,j],n=5000,level=0.95))
  }
}
ConfidenceIntervalsFemale

#In order to determine if linear models are appropriate we look at several properties of the lm
summary(lm(scaleBPnumMale[,1]~scaleBPnumMale[,3:6]))$adj.r.squared
summary(lm(scaleBPnumFemale[,1]~scaleBPnumFemale[,3:6]))$adj.r.squared


#********************************************************************************************************************************
#GLM
#For binary response variables with continuous predictors, we should use the binomial family for GLM fit
#First we create a training set and a test set for both males and females (each set is half the size of original set)
set.seed(12)
trainM<- sample(1:nrow(BPnumMale),ceiling(nrow(BPnumMale)/2))
trainF<- sample(1:nrow(BPnumFemale),ceiling(nrow(BPnumFemale)/2))
BPnumMale.train<-BPnumMale[trainM,]
BPnumMale.test<-BPnumMale[-trainM,]
BPnumFemale.train<-BPnumFemale[trainF,]
BPnumFemale.test<-BPnumFemale[-trainF,]

GLMfitM <- glm(Complications~Age+AverageBloodPressure+BiliRubinConcentration+NeutrophilPercentage,data=BPnumMale.train,family=binomial())
summary(GLMfitM) 
0.018922*mean(BPnumMale.train$Age)
0.003473*mean(BPnumMale.train$AverageBloodPressure)
0.500147*mean(BPnumMale.train$BiliRubinConcentration)
-0.273543*mean(BPnumMale.train$NeutrophilPercentage)


confint(GLMfitM) # 95% CI for the coefficients
predM<-predict(GLMfitM, BPnumMale.test, type="response") # predicted values
resM<-residuals(GLMfitM, type="deviance") # residuals 
max(predM)
min(predM)
#since the max of the predictions is less than 0.5, all the predictions will be 0 if we do not introduce a loss function
#We will classify observations as "complications = true" if the probability of complications is above 0.2
classPredM<-numeric(0)
for(i in 1:length(predM)){
  if(predM[i]>0.2) classPredM<-c(classPredM,1) else classPredM<-c(classPredM,0)
}
falsePosM<-0
truePosM<-0
trueNegM<-0
falseNegM<-0
for(i in 1:length(predM)){
  if(classPredM[i]==1 && BPnumMale.test[i,1]==0) falsePosM<-falsePosM +1
  if(classPredM[i]==1 && BPnumMale.test[i,1]==1) truePosM<-truePosM +1
  if(classPredM[i]==0 && BPnumMale.test[i,1]==1) falseNegM<-falseNegM +1
  if(classPredM[i]==0 && BPnumMale.test[i,1]==0) trueNegM<-trueNegM +1
}
falsePosM
truePosM
trueNegM
falseNegM

#in order to get the ROC and AUC we will use the pROC package
#install.packages("pROC")
library(pROC)
ROCmaleGLM <- roc(BPnumMale.test$Complications,predM,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

CImaleGLM<-ci.se(ROCmaleGLM)
plot(CImaleGLM, type="shape", col="lightblue")

#now we do the same for females:
GLMfitF <- glm(Complications~Age+AverageBloodPressure+BiliRubinConcentration+NeutrophilPercentage,data=BPnumFemale.train,family=binomial())
summary(GLMfitF) 
0.014617*mean(BPnumFemale.train$Age)
0.002936*mean(BPnumFemale.train$AverageBloodPressure)
0.411522*mean(BPnumFemale.train$BiliRubinConcentration)
-37.811879*mean(BPnumFemale.train$NeutrophilPercentage)

confint(GLMfitF) # 95% CI for the coefficients
predF<-predict(GLMfitF, BPnumFemale.test, type="response") # predicted values
resF<-residuals(GLMfitF, type="deviance") # residuals 
max(predF)
min(predF)
length(subset(predF,predF>0.5))
#since there is only 1 prediction with probability of complications above 0.5
#We will classify observations as "complications = true" if the probability of complications is above 0.2
classPredF<-numeric(0)
for(i in 1:length(predF)){
  if(predF[i]>0.2) classPredF<-c(classPredF,1) else classPredF<-c(classPredF,0)
}
falsePosF<-0
truePosF<-0
trueNegF<-0
falseNegF<-0
for(i in 1:length(predF)){
  if(classPredF[i]==1 && BPnumFemale.test[i,1]==0) falsePosF<-falsePosF +1
  if(classPredF[i]==1 && BPnumFemale.test[i,1]==1) truePosF<-truePosF +1
  if(classPredF[i]==0 && BPnumFemale.test[i,1]==1) falseNegF<-falseNegF +1
  if(classPredF[i]==0 && BPnumFemale.test[i,1]==0) trueNegF<-trueNegF +1
}
falsePosF
truePosF
trueNegF
falseNegF

#now we get the ROC curve and AUC for females:
ROCfemaleGLM <- roc(BPnumFemale.test$Complications,predF,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)

CIfemaleGLM<-ci.se(ROCfemaleGLM)
plot(CIfemaleGLM, type="shape", col="lightblue")


#**************************************************************************************************************************************
#LDA
#analysis for males
#first we seperate our data into subjects with complications and those without
numPillTrueM = BPnumMale.train[BPnumMale.train$Complications ==1,]
numPillFalseM = BPnumMale.train[BPnumMale.train$Complications ==0,]
#we remove the first 2 columns since we want the predictors only
numPill1M = numPillTrueM[,-(1:2)]
numPill0M = numPillFalseM[,-(1:2)]
#calculate the means
mu1LDAM = colMeans(numPill1M)
mu0LDAM = colMeans(numPill0M)
#compute the covariance matrix
mu1LDAmatM =  cbind( rep(mu1LDAM[1], dim(numPillTrueM)[1]), rep(mu1LDAM[2], dim(numPillTrueM)[1]), rep(mu1LDAM[3], dim(numPillTrueM)[1]), rep(mu1LDAM[4], dim(numPillTrueM)[1]) )
mu0LDAmatM =  cbind( rep(mu0LDAM[1], dim(numPillFalseM)[1]),rep(mu0LDAM[2], dim(numPillFalseM)[1]) , rep(mu0LDAM[3], dim(numPillFalseM)[1]),rep(mu0LDAM[4], dim(numPillFalseM)[1]) )
sigmaLDAM = cov(rbind(numPill1M-mu1LDAmatM, numPill0M - mu0LDAmatM))
#calculating proportion of patients who have complications:
prop1M=mean(BPnumMale.test$Complications)
prop0M=1-prop1M
#We need this library for dmvnorm function
library(mvtnorm)
#This function returns the predicted probability of having complications for a given data point:
probPredM<-function(X){
  fXcondY1LDA = dmvnorm(x=X, mean=mu1LDAM, sigma=sigmaLDAM, log=FALSE)
  fXcondY0LDA = dmvnorm(x=X, mean=mu0LDAM, sigma=sigmaLDAM, log=FALSE)
  #apply Baye's Law:
  newPred = fXcondY1LDA*prop1 / (fXcondY1LDA*prop1 + fXcondY0LDA*prop0)
  return(newPred)
}

LDAprobsM<-numeric(0)
testPredictorsM<-BPnumMale.test[,-(1:2)]
for (i in 1:length(BPnumMale.test[,1])){
  LDAprobsM<-c(LDAprobsM,probPredM(testPredictorsM[i,]))
}

ROCmaleLDA <- roc(BPnumMale.test$Complications,LDAprobsM,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)

CImaleLDA<-ci.se(ROCmaleLDA)
plot(CImaleLDA, type="shape", col="lightblue")

#now we do the LDA analysis for females:

#first we seperate our data into subjects with complications and those without
numPillTrueF = BPnumFemale.train[BPnumFemale.train$Complications ==1,]
numPillFalseF = BPnumFemale.train[BPnumFemale.train$Complications ==0,]
#we remove the first 2 columns since we want the predictors only
numPill1F = numPillTrueF[,-(1:2)]
numPill0F = numPillFalseF[,-(1:2)]
#calculate the means
mu1LDAF = colMeans(numPill1F)
mu0LDAF = colMeans(numPill0F)
#compute the covariance matrix
mu1LDAmatF =  cbind( rep(mu1LDAF[1], dim(numPillTrueF)[1]), rep(mu1LDAF[2], dim(numPillTrueF)[1]), rep(mu1LDAF[3], dim(numPillTrueF)[1]), rep(mu1LDAF[4], dim(numPillTrueF)[1]) )
mu0LDAmatF =  cbind( rep(mu0LDAF[1], dim(numPillFalseF)[1]),rep(mu0LDAF[2], dim(numPillFalseF)[1]) , rep(mu0LDAF[3], dim(numPillFalseF)[1]),rep(mu0LDAF[4], dim(numPillFalseF)[1]) )
sigmaLDAF = cov(rbind(numPill1F-mu1LDAmatF, numPill0F - mu0LDAmatF))
#calculating proportion of patients who have complications:
prop1F=mean(BPnumFemale.test$Complications)
prop0F=1-prop1F

#This function returns the predicted probability of having complications for a given data point:
probPredF<-function(X){
  fXcondY1LDA = dmvnorm(x=X, mean=mu1LDAF, sigma=sigmaLDAF, log=FALSE)
  fXcondY0LDA = dmvnorm(x=X, mean=mu0LDAF, sigma=sigmaLDAF, log=FALSE)
  #apply Baye's Law:
  newPred = fXcondY1LDA*prop1 / (fXcondY1LDA*prop1 + fXcondY0LDA*prop0)
  return(newPred)
}

LDAprobsF<-numeric(0)
testPredictorsF<-BPnumFemale.test[,-(1:2)]
for (i in 1:length(BPnumFemale.test[,1])){
  LDAprobsF<-c(LDAprobsF,probPredF(testPredictorsF[i,]))
}

ROCfemaleLDA <- roc(BPnumFemale.test$Complications,LDAprobsF,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)

CIfemaleLDA<-ci.se(ROCfemaleLDA)
plot(CIfemaleLDA, type="shape", col="lightblue")

#******************************************************************************************************************************
#QDA analysis
library(MASS)
#QDA analysis for males:
qda.M <- qda(Complications ~ Age + AverageBloodPressure + BiliRubinConcentration + NeutrophilPercentage, data = BPnumMale.train)
predQDA.M<-predict(qda.M,BPnumMale.test[,3:6])

ROCmaleQDA <- roc(BPnumMale.test$Complications,predQDA.M$posterior[,2],
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)

CImaleQDA<-ci.se(ROCmaleQDA)
plot(CImaleQDA, type="shape", col="lightblue")

#QDA analysis for females:
qda.F <- qda(Complications ~ Age + AverageBloodPressure + BiliRubinConcentration + NeutrophilPercentage, data = BPnumFemale.train)
predQDA.F<-predict(qda.F,BPnumFemale.test[,3:6])

ROCfemaleQDA <- roc(BPnumFemale.test$Complications,predQDA.F$posterior[,2],
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)

CIfemaleQDA<-ci.se(ROCfemaleQDA)
plot(CIfemaleQDA, type="shape", col="lightblue")

#*******************************************************************************************************************
#decision trees
#install.packages("tree")

#we change the logical vectors to factors so the tree function can work properly with them
BigPillData$Complications<-factor(ifelse(BigPillData$Complications,"yes","no"))
BigPillData$IsMale<-factor(ifelse(BigPillData$IsMale,"yes","no"))

#tree analysis for males
BigPillMale<-subset(BigPillData, IsMale == "yes")


BigPillMale.train<-BigPillMale[trainM,]
BigPillMale.test<-BigPillMale[-trainM,]

library(tree)
tree.BigPillDataM = tree(Complications~., data = BigPillMale,split = "gini", subset = trainM)
summary(tree.BigPillDataM)

tree.predM = predict(tree.BigPillDataM, BigPillMale.test, type = "class")
table(tree.predM,BigPillMale.test$Complications)


lossMatrix<- matrix(c(0,0.2,0.8,0),nrow=2,byrow = TRUE)

set.seed(13)
cv.BigPillMale = cv.tree(tree.BigPillDataM, FUN=prune.misclass,loss=lossMatrix)
cv.BigPillMale

prune.BigPillMale = prune.misclass(tree.BigPillDataM, best = 9)
plot(prune.BigPillMale)
text(prune.BigPillMale,pretty=0)

#tree analysis for females

BigPillFemale<-subset(BigPillData, IsMale == "no")


BigPillFemale.train<-BigPillFemale[trainF,]
BigPillFemale.test<-BigPillFemale[-trainF,]


tree.BigPillDataF = tree(Complications~., data = BigPillFemale, subset = trainF, split = "gini")
summary(tree.BigPillDataF)

tree.predF = predict(tree.BigPillDataF, BigPillFemale.test, type = "class")
table(tree.predF,BigPillFemale.test$Complications)

set.seed(13)
cv.BigPillFemale = cv.tree(tree.BigPillDataF, FUN=prune.misclass,loss=lossMatrix)
cv.BigPillFemale

prune.BigPillFemale = prune.misclass(tree.BigPillDataF, best = 9)
plot(prune.BigPillFemale)
text(prune.BigPillFemale)

#************************************************************************************************************************
#The following section is left in the code but not used in the report
#Boosting decision trees
#install.packages("gbm")
library(gbm)

#boosting analysis for males
BPmale.boost<-gbm(Complications~.-IsMale, data = BigPillMale.train, distribution = "gaussian", n.trees = 10000, shrinkage =0.01, interaction.depth = 4)
BPmale.boost
summary(BPmale.boost)

predBoostM<-predict.gbm(object = BPmale.boost,
                        newdata = BigPillMale.test,
                        n.trees = 10000,
                        type = "response")


labelsM = colnames(predBoostM)[which.max(predBoostM)]
resultM = matrix(c(BigPillMale.test$Complications, labelsM),ncol = 2, byrow = FALSE)


print(resultM)

print(pretty.gbm.tree(BPmale.boost, i.tree = BPmale.boost$n.trees))



#*********************************************************************************************************************
#The following section is left in the code but not used in the report
set.seed(12)
train = sample(1:nrow(BigPillData),2500)
BigPillData.test<-BigPillData[-train,]

library(tree)
tree.BigPillData = tree(Complications~., data = BigPillData, split = "gini", subset = train)
summary(tree.BigPillData)
tree.pred = predict(tree.BigPillData, BigPillData.test, type="class")
table(tree.pred, BigPillData.test$Complications)
(2038+34)/2500

#so our tree predicts the correct outcome 82.88% of the time

plot(tree.BigPillData)
text(tree.BigPillData)

set.seed(13)
cv.BigPillData = cv.tree(tree.BigPillData, FUN=prune.misclass)
cv.BigPillData

#so all the pruned trees misclassify 292 observations

#****************************************************************************************************************************************
#The following section is LDA for both male and female and is not used in the report, I simply did not remove it in case you wanted to have a look at it



#we select a training set half the size of the total set. Note that the seed will remain the same for the 3 methods
#to make sure the training set is the same and that we will assess the methods based on the same data
set.seed(12)
train = sample(1:nrow(BigPillData),2500)
numericBigPill.test<-numericBigPill[-train,]
numericBigPill.train<-numericBigPill[train,]
#first we seperate our data into subjects with complications and those without
numPillTrue = numericBigPill.train[numericBigPill.train$Complications ==1,]
numPillFalse = numericBigPill.train[numericBigPill.train$Complications ==0,]
#we remove the first column since we want the predictors only
numPill1 = numPillTrue[,-1]
numPill0 = numPillFalse[,-1]
#calculate the means
mu1LDA = colMeans(numPill1)
mu0LDA = colMeans(numPill0)
#compute the covariance matrix
mu1LDAmat =  cbind( rep(mu1LDA[1], dim(numPillTrue)[1]), rep(mu1LDA[2], dim(numPillTrue)[1]), rep(mu1LDA[3], dim(numPillTrue)[1]), rep(mu1LDA[4], dim(numPillTrue)[1]) , rep(mu1LDA[5], dim(numPillTrue)[1]))
mu0LDAmat =  cbind( rep(mu0LDA[1], dim(numPillFalse)[1]),rep(mu0LDA[2], dim(numPillFalse)[1]) , rep(mu0LDA[3], dim(numPillFalse)[1]),rep(mu0LDA[4], dim(numPillFalse)[1]) , rep(mu0LDA[5], dim(numPillFalse)[1]) )
sigmaLDA = cov(rbind(numPill1-mu1LDAmat, numPill0 - mu0LDAmat))

#calculating proportion of patients who have complications:
prop1=mean(numericBigPill.test$Complications)
prop0=1-prop1

#This function predicts the class of given data point and probability threshold
classPred<- function(X,prob){
  fXcondY1LDA = dmvnorm(x=X, mean=mu1LDA, sigma=sigmaLDA, log=FALSE)
  fXcondY0LDA = dmvnorm(x=X, mean=mu0LDA, sigma=sigmaLDA, log=FALSE)
  #apply Baye's Law:
  newPred = fXcondY1LDA*prop1 / (fXcondY1LDA*prop1 + fXcondY0LDA*prop0)
  if(newPred>prob){
    newPredClass = 1
  }else{
    newPredClass = 0
  }
  return(newPredClass)
}
LDApreds<-numeric(0)
testPreds<-numericBigPill.test[,-1]
for(i in 1:length(numericBigPill.test[,1])){
  LDApreds<-c(LDApreds,classPred(testPreds[i,],0.5))
}

LDALossPreds1<-numeric(0)
for(i in 1:length(numericBigPill.test[,1])){
  LDALossPreds1<-c(LDALossPreds1,classPred(testPreds[i,],0.2))
}

