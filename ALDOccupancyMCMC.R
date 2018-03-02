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


    rtruncALD=function(z,mu,tau,sigma=1){
        trunc=-mu/sigma
        if((z==0)&(trunc<=0)){
            u=runif(1)
            fn_val=-(abs(trunc)-log(u))/(1-tau)
        }else if((z==1)&(trunc>0)){
            u=runif(1)
            fn_val=trunc-log(u)/tau
        }else{
            g=pALD(q=trunc,mu=0,sigma=1,p=tau)
            if(z==1){
                min=g;max=1
            }else{
                min=0;max=g
            }
            u=min+(max-min)*runif(1)
            fn_val=qALD(prob=u,mu=0,sigma=1,p=tau)
        }
        fn_val=mu+fn_val*sigma
        return(fn_val)
    }

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
            (
                (1-pALD(phi,0,sigma,tau,TRUE))+
                (pALD(phi,0,sigma,tau,TRUE)*(1-pALD(p[,1],0,sigma,tau,TRUE))^J)
            )
        z[y0]=rbinom(sum(y0),1,round(phi.temp[y0],3))
        z1=(z==1)
        z0=(z==0)

        ##
        ## Sample v
        ##

        v=sapply(1:n,function(i)rtruncALD(
                                    z=z[i],
                                    mu=phi[i],
                                    tau=tau,
                                    sigma=sigma))

        ##
        ## Sample u
        ##

        y.mat=y[z1,]
        y.vec=c(y.mat)
        p.vec=rep(p[z1],max(J))
        u.vec=rep(NA,length(y.vec))
        u.vec=sapply(1:length(u.vec),
                     function(i)rtruncALD(z=y.vec[i],
                                          mu=p.vec[i],
                                          tau=tau,
                                          sigma=sigma))
        u[z1,]=u.vec

        ##
        ##  Sample alpha
        ##

        alpha.star=rnorm(length(alpha),alpha,alpha.tune)
        p.star=matrix(W%*%alpha.star,n,max(J))
        mh1=sum(log(dALD.m(y=u[z1,],mu=p.star[z1,],sigma=sigma,tau=tau)))
        mh2=sum(log(dALD.m(y=u[z1,],mu=p[z1,],sigma=sigma,tau=tau)))
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
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
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
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
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
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
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
        mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
        mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
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