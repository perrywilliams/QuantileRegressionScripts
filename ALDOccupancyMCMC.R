ALDOccupancyMCMC=function(data,
                          tau=0.5,
                          inits,
                          parameters=c("alpha","beta"),
                          n.iter=5000,
                          checkpoint=checkpoint,
                          output.location){
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

    alpha.mean=priors$alpha.prior[1]
    alpha.var=priors$alpha.prior[2]
    beta.mean=priors$beta.prior[1]
    beta.var=priors$beta.prior[2]

    ##
    ## Tuning parameters
    ##

    alpha.tune= 0.4237573
    beta.tune=c(0.1899601,0.2321735,0.4457457,0.1718543)
    accept.alpha=0
    accept.beta=rep(0,4)

    ##
    ## Dimensions
    ##

    n=dim(data$y)[1]
    m=length(inits$beta)
    J=data$J

    ##
    ## Pre calculations
    ##

    alpha=inits$alpha
    beta=inits$beta
    z=inits$z
    z1=(z==1)
    z0=(z==0)
    X=data$X
    W=data$W
    phi=X%*%beta
    p=W%*%alpha
    u=matrix(rnorm(n*max(J),p,1),n,max(J))
    v=matrix(rnorm(n,phi,1),n,1)
    y0=ifelse(apply(y,1,sum)==0,TRUE,FALSE)

    ##
    ## Containers
    ##

    MCMC.Chains=vector('list',length(parameters))
    names(MCMC.Chains)=parameters
    if('alpha'%in%parameters){
        MCMC.Chains$alpha=matrix(,n.iter,length(alpha))
    }
    if('beta'%in%parameters){
        MCMC.Chains$beta=matrix(,n.iter,length(beta))
    }
    if('z'%in%parameters){
        MCMC.Chains$z=matrix(,n.iter,n)
    }
    if('v'%in%parameters){
        MCMC.Chains$v=matrix(,n.iter,n)
    }
    if('phi'%in%parameters){
        MCMC.Chains$phi=matrix(,n.iter,n)
    }
    if('p'%in%parameters){
        MCMC.Chains$p=matrix(,n.iter,n)
    }

    ##
    ## Gibbs loop
    ##

    for(k in 1:n.iter){

        ##
        ## Sample z
        ##

        phi.temp=(pALD(phi,0,sigma,tau,TRUE)*
                  (1-pALD(p[,1],0,sigma,tau,TRUE))^J)/
            ((1-pALD(phi,0,sigma,tau,TRUE))+
             (pALD(phi,0,sigma,tau,TRUE)*
              (1-pALD(p[,1],0,sigma,tau,TRUE))^J)
            )
        z[y0]=rbinom(sum(y0),1,round(phi.temp[y0],3))
        z1=(z==1)
        z0=(z==0)

        ##
        ## Sample v
        ##

        v=r.truncALD(z=z,mu=phi,tau=tau,sigma=sigma)

        ##
        ## Sample u
        ##

        y.mat=y[z1,]
         y.vec=c(y.mat)
        p.vec=rep(p[z1,1],max(J))
        u.vec=rep(NA,length(y.vec))
        u.vec=r.truncALD(z=y.vec,mu=p.vec,tau=tau,sigma=sigma)
        u[z1,]=u.vec

        ##
        ##  Sample alpha
        ##

        alpha.star=rnorm(length(alpha),alpha,alpha.tune)
        p.star=matrix(W%*%alpha.star,n,max(J))
        mh1=sum(log(dALD.m(y=u[z1,],mu=p.star[z1,],sigma=sigma,tau=tau)))+
            sum(dnorm(alhpa.star,alpha.mean,sqrt(alpha.var),log=TRUE))
        mh2=sum(log(dALD.m(y=u[z1,],mu=p[z1,],sigma=sigma,tau=tau)))+
            sum(dnorm(alhpa,alpha.mean,sqrt(alpha.var),log=TRUE))

        mh=exp(mh1-mh2)
        if(mh>min(1,runif(1))){
            alpha=alpha.star
            p=p.star
            accept.alpha=accept.alpha+1
        }

        ##
        ## Sample beta0
        ##

        beta0.star=rnorm(1,beta[1],beta.tune[1])
        beta.star=c(beta0.star,beta[2:m])
        phi.star=X%*%beta.star
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))+
            dnorm(beta0.star,beta.mean,sqrt(beta.var),log=TRUE)
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))+
            dnorm(beta[1],beta.mean,sqrt(beta.var),log=TRUE)
        mh=exp(mh1-mh2)
        if(min(mh,1)>runif(1)){
            beta=beta.star
            phi=phi.star
            accept.beta[1]=accept.beta[1]+1
        }

        ##
        ## Sample beta1
        ##

        beta1.star=rnorm(1,beta[2],beta.tune[2])
        beta.star=c(beta[1],beta1.star,beta[3:m])
        phi.star=X%*%beta.star
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))+
            dnorm(beta1.star,beta.mean,sqrt(beta.var),log=TRUE)
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))+
            dnorm(beta[2],beta.mean,sqrt(beta.var),log=TRUE)
        mh=exp(mh1-mh2)
        if(min(mh,1)>runif(1)){
            beta=beta.star
            phi=phi.star
            accept.beta[2]=accept.beta[2]+1
        }

        ##
        ## Sample beta2
        ##

        beta2.star=rnorm(1,beta[3],beta.tune[3])
        beta.star=c(beta[1:2],beta2.star,beta[4:m])
        phi.star=X%*%beta.star
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))+
            dnorm(beta2.star,beta.mean,sqrt(beta.var),log=TRUE)
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))+
            dnorm(beta[3],beta.mean,sqrt(beta.var),log=TRUE)
        mh=exp(mh1-mh2)
        if(min(mh,1)>runif(1)){
            beta=beta.star
            phi=phi.star
            accept.beta[3]=accept.beta[3]+1
        }

        ##
        ## Sample beta3
        ##

        beta3.star=rnorm(1,beta[4],beta.tune[4])
        beta.star=c(beta[1:3],beta3.star)
        phi.star=X%*%beta.star
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))+
            dnorm(beta3.star,beta.mean,sqrt(beta.var),log=TRUE)
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))+
            dnorm(beta[4],beta.mean,sqrt(beta.var),log=TRUE)
        mh=exp(mh1-mh2)
        if(min(mh,1)>runif(1)){
            beta=beta.star
            phi=phi.star
            accept.beta[4]=accept.beta[4]+1
        }

        ##
        ## Save samples
        ##

        if('alpha'%in%parameters){
            MCMC.Chains$alpha[k,]=alpha
        }
        if('beta'%in%parameters){
            MCMC.Chains$beta[k,]=beta
        }
        if('z'%in%parameters){
            MCMC.Chains$z[k,]=z
        }
        if('v'%in%parameters){
            MCMC.Chains$v[k,]=v
        }
        if('phi'%in%parameters){
            MCMC.Chains$phi[k,]=phi
        }
        if('p'%in%parameters){
            MCMC.Chains$p[k,]=p[,1]
        }

        if(k%%checkpoint==0){

            ##
            ## Update tuning parameters
            ##

            if(accept.alpha/k<0.3){
                alpha.tune=alpha.tune*0.9
            }
            if(accept.alpha/k>0.5){
                alpha.tune=alpha.tune*1.1
            }

            beta.tune=ifelse(accept.beta/k>0.5,
                             beta.tune*1.1,
                      ifelse(accept.beta/k<0.3,
                             beta.tune*0.9,
                             beta.tune)
                      )

            save(MCMC.Chains,file=output.location)
        }
    }
}
