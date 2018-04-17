#########################################################################
###
###
### MCMC Algorithm to fit binomial linear mixed model using quantile
### regression
###
###
#########################################################################

ALD_Binomial_GLMM_MCMC=function(data,
                                tau=0.5,
                                inits,
                                parameters=c("alpha","beta"),
                                n.iter=5000,
                                checkpoint=1000,
                                out.loc="~/MCMC_Output.RData",
                                parallel=FALSE){
    ##
    ##  Subroutines and Packages
    ##

    library(ald)

    dALD.v=function(y,mu,sigma,tau){
        densi=ifelse(y<mu,
        (tau*(1-tau)/sigma)*exp((1-tau)*(y-mu)/sigma),
        (tau*(1-tau)/sigma)*exp(-(tau)*(y-mu)/sigma)
        )
        return(densi)
    }

    dALD.m=function(y,mu,sigma,tau){
        densi=matrix(ifelse(test=c(y)<c(mu),
        (tau*(1-tau)/sigma)*exp((1-tau)*(c(y)-c(mu))/sigma),
        (tau*(1-tau)/sigma)*exp(-(tau)*(c(y)-c(mu))/sigma)),
        dim(y)[1],dim(y)[2])
        return(densi)
    }

    p.ALD=function(q,mu=0,sigma=1,p=0.5){
        acum=ifelse(q<mu,
                    p*exp((1-p)*(q-mu)/sigma),
                    1-(1-p)*exp(-(p)*(q-mu)/sigma))
        return(acum)
    }

    q.ALD=function(prob,mu=0,sigma=1,p=0.5){
        q=ifelse(prob<p,
                 mu+((sigma*log(prob/p))/(1-p)),
                 mu-((sigma*log((1-prob)/(1-p)))/p))
        return(q)
    }

    r.truncALD=function(z,mu,tau,sigma=1){
        trunc=-mu/sigma
        g=p.ALD(q=trunc,mu=0,sigma=1,p=tau)
        min=ifelse(z==1,g,0)
        max=ifelse(z==1,1,g)
        u=min+(max-min)*runif(length(z))
        fn_val=ifelse(z==0&trunc<=0,
                      -(abs(trunc)-log(runif(length(z))))/(1-tau),
               ifelse(z==1&trunc>0,
                      trunc-log(runif(length(z)))/tau,
                      q.ALD(prob=u,mu=0,sigma=1,p=tau)
                      )
               )
        fn_val=mu+fn_val*sigma
        return(fn_val)
    }

    ##
    ## Priors
    ##

    mu.alpha.mean=priors$mu.alpha.prior[1]
    mu.alpha.var=priors$mu.alpha.prior[2]
    beta.mean=priors$beta.prior[1]
    beta.var=priors$beta.prior[2]
    var.alpha.s=priors$var.alpha.prior[1]
    var.alpha.r=priors$var.alpha.prior[2]

    ##
    ## Dimensions
    ##

    n=length(data$y)
    p=length(inits$beta)
    m=dim(data$W)[2]

    ##
    ## Data
    ##

    y=data$y
    X=data$X
    W=data$W

    ##
    ## Tuning parameters
    ##

    alpha.tune=c(0.3184822,
                 0.7391069,
                 0.9355377,
                 0.8379042,
                 0.8943194,
                 1.0344497,
                 0.9306919,
                 1.0550986,
                 0.9786558,
                 0.7731717,
                 0.9837513,
                 1.1040980)
    alpha.tune=rep(0.1,m)
    beta.tune=c(0.5,
                0.285
                )
    beta.tune=rep(0.5,p)
    accept.alpha=rep(0,m)
    accept.beta=rep(0,p)

    ##
    ## Pre calculations
    ##

    alpha=inits$alpha
    mu.alpha=inits$mu.alpha
    var.alpha=inits$var.alpha
    beta=inits$beta
    mu=X%*%beta+W%*%alpha
    v=inits$v

    ##
    ## Containers
    ##

    MCMC.Chains=vector('list',length(parameters))
    names(MCMC.Chains)=parameters
    if('alpha'%in%parameters){
        MCMC.Chains$alpha=matrix(,n.iter,m)
    }
    if('beta'%in%parameters){
        MCMC.Chains$beta=matrix(,n.iter,p)
    }
    if('v'%in%parameters){
        MCMC.Chains$v=matrix(,n.iter,n)
    }
    if('mu'%in%parameters){
        MCMC.Chains$mu=matrix(,n.iter,n)
    }
    if('mu.alpha'%in%parameters){
        MCMC.Chains$mu.alpha=matrix(,n.iter,n)
    }
    if('var.alpha'%in%parameters){
        MCMC.Chains$var.alpha=matrix(,n.iter,n)
    }
    if('tuners'%in%parameters){
        MCMC.Chains$tuners=matrix(,n.iter,m+p)
    }

    ##
    ## Gibbs loop
    ##

    for(k in 1:n.iter){

        ##
        ## Sample v
        ##

        v=r.truncALD(z=y,mu=mu,tau=tau,sigma=1)

        ##
        ##  Sample alpha
        ##

        for(i in 1:m){
            alpha.tmp=rnorm(1,alpha[i],alpha.tune[i])
            alpha.star=alpha
            alpha.star[i]=alpha.tmp
            mu.star=matrix(X%*%beta+W%*%alpha.star,n,1)
            mh1=sum(log(dALD.v(y=v,mu=mu.star,sigma=1,tau=tau)))+
                dnorm(alpha.star[i],mu.alpha,sqrt(var.alpha),log=TRUE)
            mh2=sum(log(dALD.v(y=v,mu=mu,sigma=1,tau=tau)))+
                dnorm(alpha[i],mu.alpha,sqrt(var.alpha),log=TRUE)
            mh=exp(mh1-mh2)
            if(mh>min(1,runif(1))){
                alpha=alpha.star
                mu=mu.star
                accept.alpha[i]=accept.alpha[i]+1
            }
        }

        ##
        ## Sample beta
        ##

        for(i in 1:p){
            beta.tmp=rnorm(1,beta[i],beta.tune[i])
            beta.star=beta
            beta.star[i]=beta.tmp
            mu.star=X%*%beta.star+W%*%alpha
            mh1=sum(log(dALD.v(y=v,mu=mu.star,sigma=1,tau=tau))) +
                dnorm(beta.star[i],beta.mean,sqrt(beta.var),log=TRUE)
            mh2=sum(log(dALD.v(y=v,mu=mu,sigma=1,tau=tau))) +
                dnorm(beta[i],beta.mean,sqrt(beta.var),log=TRUE)
            mh=exp(mh1-mh2)
            if(min(mh,1)>runif(1)){
                beta=beta.star
                mu=mu.star
                accept.beta[i]=accept.beta[i]+1
            }
        }

        ##
        ##  Sample mu.alpha
        ##

        A=m/var.alpha+1/mu.alpha.var
        b=(sum(alpha)/var.alpha) + (mu.alpha.mean/mu.alpha.var)
        mu.alpha=rnorm(1,mean=1/A*b,sd=sqrt(1/A))

        ##
        ##  Sample var.alpha
        ##

        a=m/2+var.alpha.s
        b=sum((alpha-mu.alpha)^2)/2 + 1/var.alpha.r
        var.alpha=1/rgamma(1,shape=a,rate=b)

        ##
        ## Save samples
        ##

        if('alpha'%in%parameters){
            MCMC.Chains$alpha[k,]=alpha
        }
        if('beta'%in%parameters){
            MCMC.Chains$beta[k,]=beta
        }
        if('v'%in%parameters){
            MCMC.Chains$v[k,]=v
        }
        if('mu'%in%parameters){
            MCMC.Chains$mu[k,]=mu
        }
        if('mu.alpha'%in%parameters){
            MCMC.Chains$mu.alpha[k,]=mu.alpha
        }
        if('var.alpha'%in%parameters){
            MCMC.Chains$var.alpha[k,]=var.alpha
        }
        if('tuners'%in%parameters){
            MCMC.Chains$tuners[k,]=c(alpha.tune,beta.tune)
        }

        ##
        ## Checkpoint
        ##

        if(k%%checkpoint==0){

            ##
            ## Update tuning parameters
            ##

            alpha.tune=ifelse(accept.alpha/k>0.5,
                              alpha.tune*1.1,
                       ifelse(accept.alpha/k<0.3,
                              alpha.tune*0.9,
                              alpha.tune)
                       )

            beta.tune=ifelse(accept.beta/k>0.5,
                             beta.tune*1.1,
                      ifelse(accept.beta/k<0.3,
                             beta.tune*0.9,
                             beta.tune)
                      )

            ##
            ## Output results
            ##


        }
    }
    return(MCMC.Chains)
}
