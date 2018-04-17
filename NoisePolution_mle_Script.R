
###
### Script to analyze the effect of noise polution on
### cavity-nesting birds
###

###
### Packages and sub-routines
###

library(lme4)
library(AICcmodavg)

###
### Load data
###

KleistData=read.csv(paste("~/Dropbox/QuantileRegression/",
               "kleist_et_al_dryad.csv",sep=""))
head(KleistData)

###
### Ash-throated fly catcher box-level analysis
###

atflpData.tmp=KleistData
atflpData=atflpData.tmp[atflpData.tmp[,6]=="ATFL",]
head(atflpData)
unoccData.tmp=KleistData
unoccData=unoccData.tmp[unoccData.tmp[,7]==0,]
atflData=rbind(atflpData,unoccData)
head(atflData)

##
## Variables
##

y=atflData$occupancy
noise=scale(atflData$noise)
fc=scale(atflData$tree_cover_box)
site=atflData$site
nb=atflData$nest_box
d=scale(atflData$distance)

##
## Models
##

m1=glmer(y~noise+fc+(1|site)+(1|nb),family="binomial")
m2=glmer(y~d+fc+(1|site)+(1|nb),family="binomial")
m3=glmer(y~fc+(1|site)+(1|nb),family="binomial")
m4=glmer(y~noise+fc+d+(1|site)+(1|nb),family="binomial")
m5=glmer(y~1+(1|site)+(1|nb),family="binomial")
summary(m1)
(aic=AIC(m1,m2,m3,m4,m5))
