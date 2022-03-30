library(MASS)
library(candisc)
require(ggplot2)
require(scales)
require(gridExtra)
library(stats)
library(foreign)
library(nnet)
library(stargazer)
library(rpart)
library(RColorBrewer)
library(rattle)

##### 625 Final Project #####
##### Investigating the relationship of body dimensions###
body<-read.table(file = "C:/Users/58257/Desktop/625 Applied multivariate analysis/Final project/body.txt")

dim(body)
colnames(body)<-c("bia","pel","bit","ched","che","elb","wri","kne","ank","shg",
                  "chg","wag","nag","hig","thg","big","fog","kng","cag","ang","wrg",
                  "age","wei","hei","gen")
body
dim(body)

bb1<-scale(body[,1:22],center = T, scale=T)
hei<-scale(body[,24],center = T,scale = T)
body5<-data.frame(bb1,hei)

body6<-data.frame(body5,body[,23])
colnames(body6)<-c("bia","pel","bit","ched","che","elb","wri","kne","ank","shg",
                        "chg","wag","nag","hig","thg","big","fog","kng","cag","ang","wrg",
                        "age","hei","wei")
body6
dim(body6)

wei<-scale(body[,23],center = T, scale = T)



body4<-data.frame(scale(body[,1:24],center = T, scale = T), body[,25])


colnames(body4)<-c("bia","pel","bit","ched","che","elb","wri","kne","ank","shg",
                   "chg","wag","nag","hig","thg","big","fog","kng","cag","ang","wrg",
                 "age","wei","hei","gen")
dim(body4)
body4
####there are 9 variables for skeletal measuremets, 12 girth measurements, 4 other variables##
### descriptive--Numerical####

### univariate analysis###


par(mfrow=c(1,1))
stars(body[1:40,1:24])


summary(body[,1:24])

apply(body[,1:24],2,var)
apply(body[,1:24],2,sd)

### let us check now for homogeneity of variance-covariance matrix assumptions.  ##

cov(body[body$gen == 1, -25])
cov(body[body$gen == 0, -25])

cor(body4[,1:24])
cov(body4[,1:24])

X <- as.matrix(body4[,c("bia","pel","bit","ched","che","elb","wri","kne","ank","shg",
                        "chg","wag","nag","hig","thg","big","fog","kng","cag","ang","wrg",
                        "wei","hei","age")])
group <- factor(body4$gen)
clinic <- lm(X ~ group)
summary(manova(clinic))

clinic1 <- manova(X ~ group)
summary(clinic1, test="Wilks")

library(car)
library(energy)
library(ICSNP)
library(MASS)

mvnorm.etest(x=body,R=199)

library(ICSNP)

HotellingsT2(body[body$gen == 1, -25], body[body$gen == 0, -25])

### principle component analysis###
library(ISLR)
library(rgl)

pca <- prcomp(X, retx=TRUE,center = TRUE, scale. = TRUE)
pca
summary(pca)
loadings <- pca$rotation
loadings 
rownames(loadings) <- colnames(X)
scores <- pca$x
scores

plot(pca)
screeplot(pca, type="lines")

dev.new(height=7, width=7)
biplot(scores[,1:8], loadings[,1:8], cex=0.7)

male<-body4[,25]==1 
plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2", col=body4[,25])
points(scores[male,1],  scores[male,2], col="blue")

library(ggbiplot)
ggbiplot(pca, choices=c(1,2), groups=body4[,25], obs.scale = 1, var.scale = 1, ellipse = T)

##### factor analysis####
### principle component solution ####
library(psych)
bb2<- principal(X, nfactors = 2, rotate = 'none', covar = TRUE, residuals=TRUE) 
bb2

bb3<- principal(X, nfactors = 8, rotate = 'none', covar = TRUE, residuals=TRUE) 
bb3

bb4<- principal(X, nfactors = 8, rotate = "varimax", covar = TRUE, residuals=TRUE) 
bb4

# principal factor solution #
bb5 <- fa(X, nfactors = 2, rotate = 'none', fm = 'pa', max.iter =1)
bb5

bb6 <- fa(X, nfactors = 8, rotate = 'none', fm = 'pa', max.iter =1)
bb6

#### MLE solution ####
bb7 <- factanal(x = X, factors = 2,  rotation = "none", method = "mle")
bb7

bb8 <- factanal(x = X, factors = 10,  rotation = "none", method = "mle")
bb8
### same result with rotation in MLE solution ####

#### multiple regressions #####


colnames(body6)<-c("bia","pel","bit","ched","che","elb","wri","kne","ank","shg",
                   "chg","wag","nag","hig","thg","big","fog","kng","cag","ang","wrg",
                   "age","hei","wei")

## regression between weight and skeletal variables & age##

mvmod <- lm(wei~ bia+pel+bit+ched+che+elb+wri+kne+ank+age+hei, data=body6)
coef(mvmod)
summary(mvmod)

anova(mvmod)
vif(mvmod)



step1 <- stepAIC(mvmod, direction="both")
step2 <- stepAIC(mvmod, direction="forward")
step3 <- stepAIC(mvmod, direction="backward")
summary(step3)

mvmod2<-lm(wei~bia+pel++bit+ched+che+wri+kne+age+hei,data=body6 )
summary(mvmod2)

## Type III SS  ##
Anova(mvmod2, type="III")

y=body6[,24]

residuals(mvmod2)

cbind(y, fitted(mvmod2), residuals(mvmod2))


## Residual Diagnostics  ###

shapiro.test(residuals(mvmod2))
qqnorm(residuals(mvmod2))

sres <- rstandard(mvmod2)
sres
sres[which(abs(sres) > 2)]

sdelres <- rstudent(mvmod2)


par(mfrow=c(1,2))
plot(sres)
plot(sdelres)

leverage <- hatvalues(mvmod2)
leverage[which(leverage > 10/32)]


cooksD <- cooks.distance(mvmod2)
influence.measures(mvmod2)


par(mfrow = c(2, 2))
plot(mvmod2)


### regression between weight and body girth variables & hei and age###
mvmod3<- lm(wei~shg+chg+wag+nag+hig+thg+big+fog+kng+cag+ang+wrg+age+hei, data=body6)
coef(mvmod3)
summary(mvmod3)
step1 <- stepAIC(mvmod3, direction="both")
step2 <- stepAIC(mvmod3, direction="forward")
step3 <- stepAIC(mvmod3, direction="backward")
summary(step3)

mvmod4<-lm(wei~shg+chg+wag+hig+thg+fog+kng+cag+age+hei,data=body6 )
summary(mvmod4)

### residual check ###

par(mfrow = c(2, 2))
plot(mvmod4)
shapiro.test(residuals(mvmod4))

#### regression all skeletal variables and girth variables plus age height on weight

mvmod5 <- lm(wei~ bia+pel+bit+ched+che+elb+wri+kne+ank+shg+chg+wag+nag+hig+thg+big+fog+kng+cag+ang+wrg+age+hei, data=body6)
summary(mvmod5)

step1 <- stepAIC(mvmod5, direction="both")
step2 <- stepAIC(mvmod5, direction="forward")
step3 <- stepAIC(mvmod5, direction="backward")
summary(step3)

mvmod6<-lm(wei~pel+ched+kne+shg+chg+wag+hig+thg+fog+kng+cag+age+hei,data=body6 )
summary(mvmod6)

### residual check ###

par(mfrow = c(2, 2))
plot(mvmod6)

#### check age variable influence on weight



library(tidyverse)

young<-filter(body, age<35)
dim(young)
young

old<-filter(body, age>36)
dim(old)
old

young<-data.frame(scale(young[,1:24]))
colnames(young)<-c( "bia","pel","bit","ched","che","elb","wri","kne","ank","shg",
                      "chg","wag","nag","hig","thg","big","fog","kng","cag","ang","wrg",
                      "age","wei","hei" )

old<-data.frame(scale(old[,1:24]))
colnames(old)<-c( "bia","pel","bit","ched","che","elb","wri","kne","ank","shg",
                         "chg","wag","nag","hig","thg","big","fog","kng","cag","ang","wrg",
                         "age","wei","hei")



mvmod7<- lm(wei ~., data=young)
summary(mvmod7)

mvmod8<- lm(wei ~ ched+kne+shg+chg+wag+hig+thg+cag+hei, data=young)
summary(mvmod8)


mvmod9<- lm(wei ~., data=old)
summary(mvmod9)

mvmod10<- lm(wei ~ ched+wag+hig+cag+age+hei, data=old)
summary(mvmod10)


### logistic regression ####

#### take first eight components input ####

body2<-scores[,1:8]

gen<-body[,25]

body3<-data.frame(body2,gen)
body3

model <- glm(gen ~ ., family = binomial, data=body3)
model

summary(model)

prob <- predict(model, body3, type="response")
round(prob, 2)

body3.pred=rep("0",507)

body3.pred[prob>0.5]<-"1"

table(body3.pred, body3[,9])

mean(body3.pred!=body3[,9])

#### drop insignificant terms
model2 <- glm(gen ~ PC1+PC2+PC3+PC5, family = binomial, data=body3)
model2

summary(model2)

prob <- predict(model2, body3, type="response")
round(prob, 2)

body3.pred=rep("0",507)

body3.pred[prob>0.5]<-"1"

table(body3.pred, body3[,9])

mean(body3.pred!=body3[,9])

#### take significant term in full regression model 

body4$gen<-as.factor(body4$gen)

model3 <- glm(gen ~pel+ched+kne+shg+chg+wag+hig+thg+fog+kng+cag+age+hei+wei , family = binomial, data=body4)

model3

summary(model3)

prob <- predict(model3, type="response")
round(prob, 2)

body.pred=rep("0",507)

body.pred[prob>0.5]<-"1"

table(body.pred, body4[,25])

mean(body.pred!=body4[,25])

#### drop terms
model4 <- glm(gen ~wag+fog+hei+wei, family = binomial, data=body4)

model4

summary(model4)

prob <- predict(model4, type="response")
round(prob, 2)

body.pred=rep("0",507)

body.pred[prob>0.5]<-"1"

table(body.pred, body4[,25])

mean(body.pred!=body4[,25])

#### classification by age groups

attach(body)


range(body[,22]) #### age range
median(body$age)

Young<-ifelse(body$age<=35,"1","2")
Young<-as.factor(Young)

newbody<-data.frame(body,Young)
newbody

model5<- glm(gen ~wag+fog+hei+wei+Young, family = binomial, data=newbody)

model5

summary(model5)

prob <- predict(model5, type="response")
round(prob, 2)

body.pred=rep("0",507)

body.pred[prob>0.5]<-"1"

table(body.pred, newbody[,25])

mean(body.pred!=newbody[,25])

## LDA  ##

library(MASS)
library(ICSNP)

body4$gen<- as.factor(body4$gen)
body.lda <- lda(gen ~ wag+fog+hei+wei, data=body4, prior=c(0.5, 0.5))
body.lda

summary(body.lda)

plot(body.lda)
predict(body.lda)$class
predict(body.lda)$posterior
table(body4$gen, predict(body.lda)$class) 
mean(body4$gen!=predict(body.lda)$class)


body.lda$means

body.lda$scaling


## cross-validated misclassification error  ##
body.lda.cv <- lda(gen ~ fog+wag+hei+wei, data=body4, CV=TRUE)
table(body4$gen, body.lda.cv$class) 
mean(body4$gen!=body.lda.cv$class)

#### Take PCA method

attach(body3)
body.lda1 <- lda(gen ~., data=body3, prior=c(0.5, 0.5))
body.lda1
summary(body.lda1)

plot(body.lda1)
predict(body.lda1)$class
predict(body.lda1)$posterior
table(body3$gen, predict(body.lda1)$class) 
mean(body3$gen!=predict(body.lda1)$class)




## QDA  ##
body.qda <- qda(gen ~ wag+fog+hei+wei, data=body4, prior=c(0.5, 0.5))
body.qda

summary(body.qda)

predict(body.qda)$class
predict(body.qda)$posterior

table(body4$gen, predict(body.qda)$class) 
mean(body4$gen!=predict(body.qda)$class)

plot(x=body4$hei, y = body4$wei, pch = as.character(body4$gen),
     col = as.character(predict(body.qda)$class))


body.qda.cv <- qda(gen ~ wag+fog+hei+wei, data=body4, CV=TRUE)
table(body4$gen, body.qda.cv$class) 
mean(body4$gen!=body.qda.cv$class)

#### Take PCA method

attach(body3)

body.qda1 <- qda(gen ~., data=body3, prior=c(0.5, 0.5))
body.qda1
table(body3$gen, predict(body.qda1)$class) 
mean(body3$gen!=predict(body.qda1)$class)


## plot LDA, QDA and PCA

ldaid <- as.integer(predict(body.lda1)$class)
qdaid <- as.integer(predict(body.qda1)$class)

par(mfrow=c(1,2))
lev <- levels(body4$gen)
plot(scores[,1:8], xlab="PC1", ylab="PC2", pch=ldaid, col=ldaid,
     main="LDA Results", xlim=c(-4, 4), ylim=c(-2, 2))

abline(h=0,lty=3)
abline(v=0,lty=3)

plot(scores[,1:8], xlab="PC1", ylab="PC2", pch=qdaid, col=qdaid,
     main="QDA Results", xlim=c(-4, 4), ylim=c(-2, 2))

abline(h=0,lty=3)
abline(v=0,lty=3)




##### K-nearest neighbors classfification #####


library(class)

set.seed(4948493) #Set the seed for reproducibility
#Sample data set (70% train, 30% test)
body.sample<-sample(1:nrow(body4),size=nrow(body4)*.7)
body.train<-body4[body.sample,] #Select the 70% of rows
body.test<-body4[-body.sample,] #Select the 30% of rows

#First Attempt to Determine Right K####
body_acc<-numeric() #Holding variable

for(i in 1:50){
  #Apply knn with k = i
  predict<-knn(body.train[,-25],body.test[,-25],
               body.train[,25],k=i)
  body_acc<-c(body_acc,
              mean(predict==body.test[,25]))
}
#Plot k= 1 through 30
plot(1-body_acc,type="l",ylab="Error Rate",
     xlab="K",main="Error Rate for body With Varying K")



#Try many Samples of Data Set to Validate K####
trial_sum<-numeric(20)
trial_n<-numeric(20)
set.seed(6033850)
for(i in 1:100){
  body_sample<-sample(1:nrow(body4),size=nrow(body4)*.7)
  body_train<-body[body_sample,]
  body_test<-body[-body_sample,]
  test_size<-nrow(body_test)
  for(j in 1:20){
    predict<-knn(body_train[,-25],body_test[,-25],
                 body_train[,25],k=j)
    trial_sum[j]<-trial_sum[j]+sum(predict==body_test[,25])
    trial_n[j]<-trial_n[j]+test_size
  }
}

plot(1-trial_sum / trial_n,type="l",ylab="Error Rate",
     xlab="K",main="Error Rate for body With Varying K (100 Samples)")



#### test wag, fog, height, weight 

body4
wag<-body4$wag
fog<-body4$fog
body7<-data.frame(wag,fog)

body8<-data.frame(body7,body4[,23])
body9<-data.frame(body8,body4[,24])
body10<-data.frame(body9,body4[,25])
body10
colnames(body10)<-c("wag","fog","wei","hei","gen")

set.seed(4948493) #Set the seed for reproducibility
#Sample data set (70% train, 30% test)
body.sample<-sample(1:nrow(body10),size=nrow(body10)*.7)
body.train<-body10[body.sample,] #Select the 70% of rows
body.test<-body10[-body.sample,] #Select the 30% of rows

### k=1 ######
body.knn<-knn(body.train[,1:4],body.test[,1:4],cl=body.train[,5],k=1)

table(body.knn,body.test[,5])
mean(body.knn!=body.test[,5])

#### k=2 ###

body.knn<-knn(body.train[,1:4],body.test[,1:4],cl=body.train[,5],k=2)
table(body.knn,body.test[,5])
mean(body.knn!=body.test[,5])

#### k=5###

body.knn<-knn(body.train[,1:4],body.test[,1:4],cl=body.train[,5],k=4)
table(body.knn,body.test[,5])
mean(body.knn!=body.test[,5])


###  cross-validation to select k  ##
knn.cv.err<-NULL
knn.cv.sd<-NULL
for (i in 1:10) { 
  temp<-NULL
  for (j in 1:1000)
    temp <- c(temp,mean(knn.cv(body.train[ ,1:4], cl = body.train$gen, k = i) != body.train$gen))
  knn.cv.err<-c(knn.cv.err,mean(temp))
  knn.cv.sd<-c(knn.cv.sd,sd(temp))
  cat("\n Done i= ",i)
}

plot(knn.cv.err, xlim = c(1, 10),  ylim=c(min(knn.cv.err - 1.96 * knn.cv.sd),
                                          max(knn.cv.err + 1.96 * knn.cv.sd)), type = "n")
lines(knn.cv.err + 1.96 * knn.cv.sd, lty = 2, col = "blue")
lines(knn.cv.err - 1.96 * knn.cv.sd, lty = 2, col = "green")
lines(knn.cv.err, col = "red")

# best choice appears to be k = 5 nearest neighbors?

body.knn.sc <- knn(body.train[,1:4],body.test[,1:4],  cl = body.train$gen, k = 5)
table(body.knn.sc,body.test[,5])
mean(body.knn.sc!=body.test[,5])






#####  CART (classification and regression tree)
library(tree)

body4$gen<-as.factor(body4$gen)

body.tree1 <- tree(gen ~., data = body4)

body.tree1
summary(body.tree1)

plot(body.tree1)
text(body.tree1)

# prune the tree  #
body.tree.cv <- cv.tree(body.tree1, K = nrow(body4))
plot(body.tree.cv)

body.prune.tree<-prune.tree(body.tree1, best = 8)

plot(body.prune.tree)
text(body.prune.tree)

summary(body.prune.tree)

## different split rule  ##
body.tree2 <- tree(gen ~ ., data = body4, split = "gini")
body.tree2
summary(body.tree2)

plot(body.tree2)
text(body.tree2)


# colorful tree  #

library(rpart)
library(ggplot2)
library(RColorBrewer)
library(rattle)

body.rpart <- rpart(gen ~ ., data=body4, method="class")

body.rpart
summary(body.rpart)

fancyRpartPlot(body.rpart, main="body build")


##### Tree in sample by seperating age into two groups

attach(body)


range(body[,22]) #### age range
median(body$age)

Young<-ifelse(body$age<=35,"1","2")
Young<-as.factor(Young)

newbody<-data.frame(body,Young)
newbody

library(tree)
body.tree3 <- tree(Young ~.-age, data = newbody)

body.tree3
summary(body.tree3)

plot(body.tree3)
text(body.tree3)

# prune the tree  #
body.tree.cv <- cv.tree(body.tree3, K = nrow(newbody))
plot(body.tree.cv)

body.prune.tree<-prune.tree(body.tree3, best = 3)

plot(body.prune.tree)
text(body.prune.tree)

summary(body.prune.tree)

## different split rule  ##
body.tree4 <- tree(Young ~. -age, data = newbody, split = "gini")
body.tree4
summary(body.tree4)

plot(body.tree4)
text(body.tree4)


# colorful tree  #

library(rpart)
library(ggplot2)
library(RColorBrewer)
library(rattle)

body.rpart <- rpart(Young ~. -age, data=newbody, method="class")

body.rpart
summary(body.rpart)

fancyRpartPlot(body.rpart, main="body build")

#### Take PCA method ####

attach(body3)


body.tree5 <- tree(gen ~., data=body3)
body.tree5
summary(body.tree5)

plot(body.tree5)
text(body.tree5)

predict(body.tree5, newdata=body3)






##############################################################
######  Bagging
##############################################################

library(ipred)


body.bag1 <- bagging(as.factor(gen) ~ ., data = body4, coob = T)
body.bag1 
summary(body.bag1)

predict(body.bag1)

acc=(body[,25]==predict(body.bag1))
1-length(acc[acc=="TRUE"])/length(acc)

body.bag1$mtrees[[1]]

body.bag2 <- bagging(as.factor(gen) ~ ., data = body4, nbagg=100, coob = T)
body.bag2 
summary(body.bag2)

predict(body.bag2)

acc=(body[,25]==predict(body.bag2))
1-length(acc[acc=="TRUE"])/length(acc)




library(rpart)
body.bag3 <- bagging(as.factor(gen) ~ ., data = body4, 
                     control=rpart.control(maxdepth=5),nbagg=100, coob = T)
body.bag3
summary(body.bag3)
body.bag3$mtrees[[1]]

###### do bagging in age two groups

library(ipred)


body.bag4 <- bagging(as.factor(Young) ~.-age, data = newbody, coob = T)
body.bag4

predict(body.bag4)

acc=(newbody[,26]==predict(body.bag4))
1-length(acc[acc=="TRUE"])/length(acc)

body.bag4$mtrees[[1]]

body.bag5 <- bagging(as.factor(Young) ~.-age, data = newbody, nbagg=100, coob = T)
body.bag5


predict(body.bag2)

acc=(body[,25]==predict(body.bag2))
1-length(acc[acc=="TRUE"])/length(acc)




library(rpart)
body.bag3 <- bagging(as.factor(gen) ~ ., data = body4, 
                     control=rpart.control(maxdepth=5),nbagg=100, coob = T)
body.bag3
summary(body.bag3)
body.bag3$mtrees[[1]]

#### take PCA method

body.bag6 <- bagging(as.factor(gen) ~., data = body3, coob = T)
body.bag6
summary(body.bag6)



##############################################################
######  Boosting
##############################################################

library(gbm)

body.boost1 <- gbm(gen ~ ., data = body4, distribution = "adaboost", n.trees = 100)
summary(body.boost1)

par(mfrow=c(1,2))
plot(body.boost1,i="fog")
plot(body.boost1,i="shg")



best.iter <- gbm.perf(body.boost1, method="OOB")
print(best.iter)
summary(body.boost1, n.trees=best.iter) 

confusion <- function(a, b){
  tbl <- table(a, b)
  mis <- 1 - sum(diag(tbl))/sum(tbl)
  list(table = tbl, misclass.prob = mis)
}

body.predict <- (predict(body.boost1, body4, n.trees=best.iter)>0)*TRUE 
confusion(body.predict, body4$gen)


body.boost2 <- gbm(gen ~ ., data = body4, distribution = "adaboost", n.trees = 500)
summary(body.boost2)

best.iter <- gbm.perf(body.boost2, method="OOB")
print(best.iter)
summary(body.boost2, n.trees=best.iter) 

body.predict <- (predict(body.boost2, body4, n.trees=best.iter)>0)*TRUE 
confusion(body.predict, body4$gen)

body.predict <- (predict(body.boost2, body4, n.trees=400)>0)*TRUE  # boosting resistant to overfitting
confusion(body.predict, body4$gen)


### take PCA method

body.boost3 <- gbm(gen ~ ., data = body3, distribution = "adaboost", n.trees = 100)
summary(body.boost3)

body.predict <- (predict(body.boost3, body3, n.trees=best.iter)>0)*TRUE 
confusion(body.predict, body3$gen)


body.boost4 <- gbm(gen ~ ., data = body3, distribution = "adaboost", n.trees = 500)
summary(body.boost4)

body.predict <- (predict(body.boost4, body3, n.trees=best.iter)>0)*TRUE 
confusion(body.predict, body3$gen)


##############################################################
#   Random forests
##############################################################

library(randomForest)

body.rf <- randomForest(gen ~ ., data = body4)

confusion(predict(body.rf, body4), body4$gen)
body.rf$err.rate

plot(body.rf)

getTree(body.rf, 1, labelVar=TRUE)

varImpPlot(body.rf)

predict(body.rf,type="prob")

#  Importance
importance(body.rf, type=2)



##
body.rf1 <- randomForest(gen ~ ., data = body4, mtry=8, importance=TRUE, do.trace=400)
body.rf1



#  Try on all variables over 10
body.rf2 <- randomForest(gen ~ bia+elb+shg+chg+big+fog+wrg, data = body4) 

confusion(predict(body.rf2, body4), body4$gen)


#  Try on top 5
body.rf3 <- randomForest(gen ~ fog+bia+shg+wrg+elb, data = body4) 

confusion(predict(body.rf3, body4), body4$gen)

##################################
#  trace plot of OOB error estimate as function of size of forest

plot(body.rf2, log="y")

#  Importance plot:
varImpPlot(body.rf2)

par(mfrow=c(1,2))
varImpPlot(body.rf2, class="0")
varImpPlot(body.rf2, class="1")

## the effect of different values of a given variable on the class prediction
par(mfrow=c(1,2))
partialPlot(body.rf2, body4, fog, "1")
partialPlot(body.rf2, body4, fog, "0")

### Take PCA method
body.rf4 <- randomForest(gen ~ ., data = body3, mtry=8, importance=TRUE, do.trace=100)
body.rf4

plot(body.rf4)

getTree(body.rf4, 1, labelVar=TRUE)

varImpPlot(body.rf4)

#  Importance
importance(body.rf4, type=2)



