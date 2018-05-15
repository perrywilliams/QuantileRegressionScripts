#########################################################################
###
### ALD Occupancy Model
###
### Created 01Mar2018
###
#########################################################################

rm(list=ls())

###
### Packages and sub-routines
###

required.packages=c("ald",
                    "bayesQR",
                    "ggplot2",
                    "RCurl",
                    "parallel"
                    )
lapply(required.packages,library,character.only=TRUE)

## Wrapper for parallel processing
run.chain.2pl.list=function(tau=tau.v,
                            data=data,
                            inits=inits,
                            parameters=parameters,
                            n.iter=n.iter,
                            checkpoint=checkpoint,
                            out.loc=out.loc.l
                            ){
    chain.list=mclapply(1:length(tau),
                        function(j){
                            ## Set variables for this core
                            this.tau=tau[j]
                            this.out.loc=out.loc[[j]]

                            ## Run the chain on this core
                            this.chain=
                                ALD_Binomial_GLM_MCMC(data=data,
                                                       tau=this.tau,
                                                       inits=inits,
                                                       parameters=
                                                           parameters,
                                                       n.iter=n.iter,
                                                       checkpoint=
                                                           checkpoint,
                                                       out.loc=
                                                           this.out.loc
                                                       )
                            return(this.chain)
                        },
                        mc.cores=min(length(tau.v),detectCores())
                        )

    ## Save the initial random seed as the name of the chain
    names(chain.list)=tau.v
    return(chain.list)
}

###
### calc CDF of ALD
###

p.ALD=function(q,mu=0,sigma=1,p=0.5){
        acum=ifelse(q<mu,
                    p*exp((1-p)*(q-mu)/sigma),
                    1-(1-p)*exp(-(p)*(q-mu)/sigma))
        return(acum)
}

###
### Quantiles
###

tau.v=c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9)
## tau.v=0.5

###
### Load real data
###

ThreshData=read.csv(paste0("/Users/pwill/Dropbox/Projects/",
                "QuantileRegression/3552231/",
                "Supplement_Avian_data.csv"))
ThreshData=ThreshData[!is.na(ThreshData$AGE),]

###
### Species:  MGWA (no), OCWA (yes:broad,age), SWTH (yes:con,age),
###           WIFL (no), WIWA (no)
###

species="OCWA"

###
### Data
###

y=apply(ThreshData[,grepl(species,names(ThreshData))],1,max)
y=matrix(y,,1)
J=apply(!is.na(y),1,sum)
n=dim(y)[1]

###
### Occupancy Covariates
###

con=(ThreshData$CONIFER)
broad=(ThreshData$BROAD)
dbroad=(ThreshData$DECBROAD)
hwd=(ThreshData$HWD)
hwd2=(ThreshData$HWD2000M)
ele=(ThreshData$ELEVM)
age=(ThreshData$AGE)

###
### Occupancy Model
###

X=cbind(con, broad, dbroad, hwd, ele, age)     #threshold model
m=dim(X)[2]

###
### Fit with bayesQR()
###

quantvec=seq(0.1,0.9,0.1)
out=bayesQR(y~X,quantile=quantvec,alasso=TRUE,ndraw=30000)
summary(out,burnin=5000)
plot(out,plottype="trace",burnin=15000)

plot(out,plottype="quantile",burnin=15000)



ind=1:9
bar=summary(out[[1]])[[1]]$betadraw[2,]
for(i in 2:9){
    out[[i]]=bayesQR(y~X,quantile=ind[i]/10,ndraw=n.iter)
    foo=summary(out[[i]])[[1]]$betadraw[2,]
    bar=rbind(bar,foo)
}

plot(seq(.1,.9,.1),bar[,1],pch=16,ylim=c(-20.5,20.5))
segments(seq(.1,.9,.1),bar[,2],seq(.1,.9,.1),bar[,1])
segments(seq(.1,.9,.1),bar[,3],seq(.1,.9,.1),bar[,1])
abline(h=0)

###
### Detection Covariates
###

W=matrix(1,n,1)

###
### Detection Model
###

W=W
q=dim(W)[2]

## ######################################################################
## ### Alternatively, simulate occupancy data
## ######################################################################

## n=500
## m=4
## q=4
## alpha.truth=rnorm(m)
## beta.truth=rnorm(q)
## sigma=1
## J=rep(3,n)
## x1=rnorm(n)
## x2=rnorm(n)
## x3=rnorm(n)
## X=cbind(1,x1,x2,x3)
## w1=rnorm(n)
## w2=rnorm(n)
## w3=x2^2
## W=cbind(1,w1,w2,w3)
## psi.truth=1-p.ALD(-X%*%beta.truth,mu=0,sigma=1,p=0.5)
## z.truth=rbinom(n,1,psi.truth)
## p.truth=1-p.ALD(-W%*%alpha.truth,mu=0,sigma=1,p=0.5)
## y=matrix(0,n,max(J))
## for(i in 1:max(J)){
##      y[z.truth==1,i]=rbinom(sum(z.truth),1,p.truth[z.truth==1])
## }

#########################################################################

###
### Bundle data for MCMC sampler
###

data=list(y=y,
          X=X## ,
          ## W=W,
          ## J=J
          )

###
### MCMC Settings
###

n.iter=30000
checkpoint=100

###
### Prior hyperparameters
###

## alpha
alpha.mean=0
alpha.var=10

## beta
beta.mean=0
beta.var=10

priors=list(alpha.prior=c(alpha.mean,alpha.var),
            beta.prior=c(beta.mean,beta.var)
            )

###
### Starting Values
###

alpha=matrix(rnorm(q),,1)
beta=matrix(rnorm(m),,1)
z=matrix(ifelse(apply(y,1,sum)>0,1,0),n,1)

inits=list(alpha=alpha,
           beta=beta,
           z=z
           )

###
### Parameters to monitor
###

parameters=c("alpha",
             "beta",
             "tuners",
             "t")

###
### Output location
###

out.loc.l=list()
for(i in 1:length(tau.v)){
    out.loc.l[[i]]=paste("~/Dropbox/Projects/QuantileRegression/",
                         "Threshold/",
                         species,
                         "/MCMC.Output.",
                         tau.v[i],
                         ".RData",sep="")
}


###
### Source 'ALDOccupancyMCMC' function (locally or GitHub)
###

## script=getURL(
##     paste("https://raw.githubusercontent.com/",
##           "perrywilliams/QuantileRegressionScripts/",
##           "master/ALDOccupancyMCMC2.R",sep=""),
##     ssl.verifypeer=FALSE)
## eval(parse(text = script))

source(paste("~/Dropbox/GitHub/",
             "QuantileRegression-master/",
             "QuantileRegression/ALD_Binomial_GLM_Threshold.R",
             sep=""))

#########################################################################
### Fit model with 'ALDOccupancyMCMC' function
#########################################################################

sys.time=Sys.time()
MCMCOutput=run.chain.2pl.list(tau=tau.v,
                              data=data,
                              inits=inits,
                              parameters=parameters,
                              n.iter=n.iter,
                              checkpoint=checkpoint,
                              out.loc=out.loc.l
                              )
Sys.time()-sys.time
save(MCMCOutput,file=paste("~/Dropbox/Projects/QuantileRegression/",
                           "ModelOutput/Threshold/",
                           species,
                           "/MCMC.Output.RData",
                           sep=""))

#########################################################################
### Load output
#########################################################################

load(paste("~/Dropbox/Projects/QuantileRegression/",
           "ModelOutput/Threshold/",
           species,
           "/MCMC.Output.RData",
           sep=""))

(status=sum(!is.na(MCMCOutput[[1]]$beta[,1])))
(burn=round(status/10))
thin=1
ind=seq(burn,status,thin)
length(ind)

#########################################################################
### Check convergence
#########################################################################

###
### Trace plots of alpha
###

par(mfrow=c(6,2),mar=c(4,4,1,1))
for(j in 1:length(tau.v)){
    for(i in 1:m){
        plot(MCMCOutput[[j]]$alpha[ind,i],
             type='l',
             ylab="")
         ## abline(h=alpha.truth[i],col=2)
        plot(MCMCOutput[[j]]$tuner[ind,i],
             type='l',
             ylab="")
        print(i)
    }
}

###
### Trace plots of beta
###

par(mfrow=c(3,2),mar=c(4,4,1,1))
for(j in 1:length(tau.v)){
    for(i in 1:m){
        plot(MCMCOutput[[j]]$beta[,i],type='l',
             ylab="")
        ## abline(h=beta.truth[i],col=2)
        plot(MCMCOutput[[j]]$tuner[,i],
             type='l',
             ylab="")
    }
    readline()
}


###
### Trace plots of t
###

plot(MCMCOutput[[1]]$t,type='l')


#########################################################################
### Quantile regression inference
#########################################################################

beta=1

betas=matrix(,length(ind),length(tau.v))
for(i in 1:length(tau.v)){
    betas[,i]=MCMCOutput[[i]]$beta[ind,beta]
}

###
### Violin plots
###

foo=cbind(betas[,1],1/10)
for(i in 2:length(tau.v)){
    bar=cbind(betas[,i],tau.v[i])
    foo=rbind(foo,bar)
}
bar=as.data.frame(foo)
names(bar)=c("beta","tau")
bar$tau=as.factor(bar$tau)
ggplot(bar, aes(x=tau, y=beta)) +
    geom_violin()

###
### Look at prob(occ|tau)
###

rm(foo)
tau=0.5

beta.tau=MCMCOutput[[tau*10]]$beta[ind,]
X=cbind(1,0)
mu=X%*%t(beta.tau)
p.occ=1-p.ALD(q=-mu,mu=0,sigma=1,p=tau)
i=0
foo=data.frame(i,t(p.occ),tau)
for(i in 2:18){
    X=cbind(1,i)
    mu=X%*%t(beta.tau)
    p.occ=1-p.ALD(q=-mu,mu=0,sigma=1,p=tau)
    bar=data.frame(i,t(p.occ),tau)
    foo=rbind(foo,bar)
}
foo=as.data.frame(foo)
names(foo)=c("X","Y","tau")
foo$X=as.factor(foo$X)

ggplot(foo, aes(x=X, y=Y)) +
    geom_violin()+
    ylim(0,1)


b.m=apply(betas,2,mean)
b.l=apply(betas,2,quantile,0.025)
b.u=apply(betas,2,quantile,.975)

x.axs=c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9)
plot(x.axs,b.m,
     ylim=c(-5,5))
abline(h=0)
segments(x.axs,b.l,x.axs,b.m)
segments(x.axs,b.u,x.axs,b.m)
