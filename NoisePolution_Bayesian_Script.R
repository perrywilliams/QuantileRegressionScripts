#############################################################################
###
### Script to analyze the effect of noise polution on
### cavity-nesting birds using binomial linear mixed model and
### Bayesian quantile regression
### Uses function "ALD_Binomial_GLMM_MCMC.R"
###
#############################################################################

rm(list=ls())

###
### Packages and sub-routines
###

required.packages=c("ald",
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
                            out.loc=out.loc.l,
                            parallel=TRUE
                            ){
    chain.list=mclapply(1:length(tau),
                        function(j){
                            ## Set variables for this core
                            this.tau=tau[j]
                            this.out.loc=out.loc[[j]]

                            ## Run the chain on this core
                            this.chain=
                                ALD_Binomial_GLMM_MCMC(data=data,
                                                       tau=this.tau,
                                                       inits=inits,
                                                       parameters=
                                                           parameters,
                                                       n.iter=n.iter,
                                                       checkpoint=
                                                           checkpoint,
                                                       out.loc=
                                                           this.out.loc,
                                                       parallel=
                                                           parallel
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
### Quantiles
###

tau.v=c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9)
## tau.v=0.5

###
### Species ("JUTI","BEWR","ATFL","WEBL","MOBL")
###

## species="MOBL" done (200k iterations)
species="ATFL"

###
### Load real data
###

KleistData=read.csv(paste("~/Dropbox/QuantileRegression/",
               "kleist_et_al_dryad.csv",sep=""))

###
### Data
###

SpeciesData.tmp=KleistData
SpeciesData=SpeciesData.tmp[SpeciesData.tmp[,6]==species,]
unoccData.tmp=KleistData
unoccData=unoccData.tmp[unoccData.tmp[,7]==0,]
SpeciesData=rbind(SpeciesData,unoccData)

###
### Response Variable
###

y=SpeciesData$occupancy
n=length(y)

###
### 'Fixed' effects
###

noise=scale(SpeciesData$noise)
fc=scale(SpeciesData$tree_cover_box)
d=scale(SpeciesData$distance)

## Model
X=cbind(noise,fc)

###
### 'Random' effects
###

site=SpeciesData$site
nb=SpeciesData$nest_box
pair=as.factor(SpeciesData$pair)
bs=as.factor(SpeciesData$box_status)

## Model
W=model.matrix(y~pair)

###
### Dimensions
###

p=dim(X)[2]
m=dim(W)[2]

#############################################################################
### Simulate data with known parameters to evaluate model
#############################################################################

set.seed(2018)
mu.alpha.truth=5
var.alpha.truth=1.5^2
alpha.truth=rnorm(m,mu.alpha.truth,sqrt(var.alpha.truth))
beta.truth=rnorm(p,0,sqrt(10^2))
mu.truth=X%*%beta.truth+W%*%alpha.truth
v.truth=sapply(1:n,function(i)rALD(1,mu=mu.truth[i],sigma=1,p=tau.v))
y=ifelse(v.truth>0,1,0)

#############################################################################

###
### Bundle data for MCMC sampler
###

data=list(y=y,
          X=X,
          W=W
          )

###
### MCMC Settings
###

n.iter=100000
checkpoint=100

###
### Prior hyperparameters
###

## beta
beta.mean=0
beta.var=100^2

## mu.alpha
mu.alpha.mean=0
mu.alpha.var=10^2

## var.alpha
var.alpha.s=1
var.alpha.r=0.5

priors=list(mu.alpha.prior=c(mu.alpha.mean,mu.alpha.var),
            beta.prior=c(beta.mean,beta.var),
            var.alpha.prior=c(var.alpha.s,var.alpha.r)
            )

###
### Starting Values
###

beta=c(-0.4755743, -0.5875835)
beta=rnorm(dim(X)[2],beta.mean,1)
mu.alpha=rnorm(1,mu.alpha.mean,sqrt(mu.alpha.var))
var.alpha=1/rgamma(1,shape=var.alpha.s,rate=var.alpha.r)
alpha=c(-17.242589,
        -10.036379,
        -1.549096,
        -19.083194,
        -1.950155,
        -9.867848,
        -11.497288,
        -7.766612,
        -6.752701,
        -12.906468,
        3.137763,
        -4.679458)
alpha=rnorm(dim(W)[2],mu.alpha,sqrt(var.alpha))

mu=X%*%beta+W%*%alpha
v=sapply(1:n,function(i)rALD(1,mu=mu[i],sigma=1,p=tau.v))

inits=list(alpha=alpha,
           beta=beta,
           v=v,
           mu.alpha=mu.alpha,
           var.alpha=var.alpha
           )

###
### Parameters to monitor
###

parameters=c("alpha",
             "beta",
             "mu.alpha",
             "var.alpha",
             "tuners"
             )

###
### Output location
###

species="Simulation"
out.loc.l=list()
for(i in 1:length(tau.v)){
    out.loc.l[[i]]=paste("~/Dropbox/QuantileRegression/",
                         "ModelOutput/Box/",
                         species,
                         "/MCMC.Output.",
                         tau.v[i],
                         ".RData",sep="")
}

###
### Source 'ALD_Binomial_GLMM_MCMC' function (locally or GitHub)
###

## script=getURL(
##     paste("https://raw.githubusercontent.com/",
##           "perrywilliams/QuantileRegressionScripts/",
##           "master/ALDOccupancyMCMC.R",sep=""),
##     ssl.verifypeer=FALSE)
## eval(parse(text = script))

source(paste("~/Dropbox/GitHub/",
             "QuantileRegression-master/",
             "QuantileRegression/ALD_Binomial_GLMM_MCMC.R",
             sep=""))

#############################################################################
### Fit model with 'ALDOccupancyMCMC' function
#############################################################################

sys.time=Sys.time()
MCMCOutput=run.chain.2pl.list(tau=tau.v,
                              data=data,
                              inits=inits,
                              parameters=parameters,
                              n.iter=n.iter,
                              checkpoint=checkpoint,
                              out.loc=out.loc.l,
                              parallel=FALSE
                              )
Sys.time()-sys.time
save(MCMCOutput,file=paste("~/Dropbox/QuantileRegression/",
                           "ModelOutput/Box/",
                           species,
                           "/MCMC.Output.RData",
                           sep=""))

#############################################################################
### Load output
#############################################################################

species="ATFL"
load(paste("~/Dropbox/QuantileRegression/",
           "ModelOutput/Box/",
           species,
           "/MCMC.Output.RData",
           sep=""))

(status=sum(!is.na(MCMCOutput[[1]]$alpha[,1])))
(burn=round(status/10))
thin=10
ind=seq(burn,status,thin)
length(ind)
acf(MCMCOutput[[1]]$beta[ind,1])

#############################################################################
### Check convergence
#############################################################################

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
   ## readline()
}

###
### Trace plots of beta
###

par(mfrow=c(2,2),mar=c(4,4,1,1))
for(j in 1:length(tau.v)){
    for(i in 1:p){
        plot(MCMCOutput[[j]]$beta[,i],type='l',
             ylab="")
        ## abline(h=beta.truth[i],col=2)
        ## plot(MCMCOutput[[j]]$tuner[,m+i],
        ##      type='l',
        ##      ylab="")
    }
    readline()
}

###
### Trace plots of mu.alpha
###

par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(MCMCOutput[[1]]$mu.alpha[ind,1],type='l',
     ylab="")

###
### Trace plots of var.alpha
###

plot(MCMCOutput[[1]]$var.alpha[ind,1],type='l',
     ylab="")

#############################################################################
### Quantile regression inference
#############################################################################

beta=2
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

