rm(list=ls())
library(bayesQR)
library(ald)

Data=read.csv("~/Dropbox/Projects/QuantileRegression/wtmatrix.csv")
head(Data)
y=Data[,1]## apply(Data[,1:3],1,max,na.rm=TRUE) ##Data[,2]##
X=cbind(scale(Data$forest),scale(Data$elev),scale(Data$elev)^2)
q.v=seq(0.05,0.95,0.05)
out.wt=bayesQR(y~X,quantile=q.v,alasso=FALSE,ndraw=15000)

## plot(out.wt,plottype="trace")
plot(out.wt,plottype="quantile")

X.tmp=seq(range(X)[1],range(X)[2],0.1)
X.predict=cbind(1,
                0,
                X.tmp,
                X.tmp^2
                )
foo=predict(out.wt,X.predict)
plot(seq(range(X)[1],range(X)[2],0.1),
     foo,type='l')

## Calculate predicted probabilities
pred=predict(object=out.wt, X=cbind(1,X), burnin=5000)

## Make histogram of predicted probabilities
hist(pred,breaks=10)

## Calculate Percentage Correclty Classified (PCC)
mean(y==as.numeric(pred>.5))

## Compare with logit model
mylogit <- glm(y ~ X, family=binomial(logit))

## Make histogram of predicted probabilities
hist(mylogit$fit,breaks=10)

## Calculate Percentage Correclty Classified (PCC)
mean(y==as.numeric(mylogit$fit>.5))

###
###
###

coefs=summary(out.wt)[[1]]$betadraw[,1]
plot(X.tmp,
     1-pALD(-(coefs[1]+coefs[3]*X.tmp+coefs[4]*X.tmp^2),
            mu=0,p=1/20),type='l',
     ylim=c(0,0.65),ylab="P(y=1|X)",
     xlab="Elevation"
     )

for(i in 2:length(summary(out.wt))){
    coefs=summary(out.wt)[[i]]$betadraw[,1]
    lines(X.tmp,
          1-pALD(-(coefs[1]+coefs[3]*X.tmp+coefs[4]*X.tmp^2),mu=0,p=i/20),
          type='l',col=i)
    readline()
}
points(X[,2],y*.64,pch=16)
