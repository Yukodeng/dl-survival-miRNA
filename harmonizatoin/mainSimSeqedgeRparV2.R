#### Simulate data with batch effects based on the edgeR method and simulated expression data #####
#### Batch effects on both mean and realized counts ####
#### Batch effect variance is proportional to gene expression variance so batch effect is different for different genes ####

mainSimSeqedgeRpar=function(beta0, nonzero_position, n_train, n_test, ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size, 
                             he_train, he_test, batch_sd_mu, batch_sd_obs, batch_sort_train, batch_sort_test, stratify_train, stratify_test, cmbt_train, cmbt_test, cmbt_frozen, norm_type, 
                             norm_train, norm_test, standardize=0, univ_cutoff_range=c(0,0.05), nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=6, 
                             addon, nsim, seed){

  library(survival)
  library(preprocessCore)
  library(glmnet)
  library(sva)
  library(MASS)
  library(stringr)
  library(edgeR)
  library(R.utils)

  options(stringsAsFactors=F)
  setwd("/users/PAS1476/andyni/Augseq/")
  # setwd("C:/Users/ni.304/OneDrive - The Ohio State University/Research/R21 Nov2020 submission with Li-Xuan Qin/Augmented sequencing data analysis/")
  source("norm.R")
  source("glmnet_strat_cv.R")
  
  cat(paste(addon,"-seed-",seed))
  
  betas.gene=read.csv("betas.gene.462.csv")$x  ## estimated gene-specific intercept
  libsizes=read.csv("libsizes.462.csv")$x  ## estimated library size
  disps.gene=read.csv("disps.gene.462.csv")$x  ## estimated gene-specific dispersion factor
  corr=read.csv("cor.mat.462.csv")  ## estimated correlation matrix among genes
  geneid=names(corr)
  corr=as.matrix(corr)
  p=length(geneid)
  
  set.seed(seed)
  
  # letter after "_": o - oracle analysis, u - univariate analysis, l - lasso-penalized reg, a - adaptive lasso-penalized reg
  
  b_o=b_u=b_l=b_a=matrix(NA, nrow=nsim, ncol=p)
  c_o_test=c_u_test=c_l_test=c_a_test=l_l=l_a=rep(NA, nsim)
  subsize=rep(NA, nsim)

  if(stratify_train==1){
    b_os=b_us=b_ls=b_as=matrix(NA, nrow=nsim, ncol=p)
    c_os_test=c_us_test=c_ls_test=c_as_test=rep(NA, nsim)
    l_ls=l_as=al_ls=al_as=rep(NA, nsim)
  }
  
  if(cmbt_train==1){
    b_oc=b_uc=b_lc=b_ac=matrix(NA, nrow=nsim, ncol=p)
    c_oc_test=c_uc_test=c_lc_test=c_ac_test=rep(NA, nsim)
    l_lc=l_ac=al_lc=al_ac=rep(NA, nsim)
    subsize_c=rep(NA, nsim)
  }
  
  for(iter in 1:nsim){
    tryCatch(expr = {withTimeout({
      
    ## training data
    # generate gamma-distributed true gene expression 
    norm.train=mvrnorm(n_train, mu=rep(0,p), Sigma=corr)
    unif.train=pnorm(norm.train)   
    gamma.train=matrix(NA, nrow=n_train, ncol=p)
    for(g in 1:p){
      gamma.train[,g]=qgamma(unif.train[,g], shape=1/disps.gene[g], scale=exp(betas.gene[g])*disps.gene[g])
    }
    
    # generate survival time based on gamma-distributed true gene expression
    xtrue.train=gamma.train*2*10^6  ## scale to observed counts
    lxtrue.train=log2(xtrue.train+0.5)
    sds.lxtrue.train=apply(lxtrue.train, 2, sd)
    # h0=1.5  ## for small gene effect size
    h0=2200  ## for large gene effect size
    b0=rep(0,p)
    b0[nonzero_position]=beta0
    h=h0*exp(lxtrue.train%*%b0)
    t0=rexp(n_train, h)
    ct=runif(n_train, 0, 300) 
    t=ifelse(t0<=ct, t0, ct)
    delta=1*(t0<=ct)  ## about 75% event rate
    
    id=1:n_train
    ncase.train=sum(delta)

    # sample library size (i.e. sequencing depth) from the observed library sizes of the original MSK sarcoma dataset
    temp.sample=sample(libsizes, n_train, replace=T)
    libsizes.sample=temp.sample+round(runif(n_train, -10000, 10000))  # add some perturbation to avoid identical library sizes
    libsizes.sample.vec=rep(libsizes.sample, p)   
    # libsizes.sample.vec=rep(10^6, n_train*p) 
    
    # introduce batch effect
    batch.id=rep(1:(n_train/batch_size+1), each=batch_size)[1:n_train]
    batch.id.unique=unique(batch.id)
    if(he_train==1){
      batch.n=length(batch.id.unique)
      # simulate batch effects on the scale parameter of gamma distribution (i.e. batch effect on  mean and variance of negative binomial counts)
      batch.mu=mvrnorm(batch.n, rep(0,p), (batch_sd_mu*diag(sds.lxtrue.train))^2)  ## batch sd is proportional to gene sd
      batch.mu=apply(batch.mu, 2, sort, decreasing=T)  ## sort batch effect for each gene
      batch.mu.mat=batch.mu[rep(1:batch.n, each=batch_size),]
      # simulate batch effects on the realized counts from negative binomial distribution, which correlate with batch effects on gamma distribution
      batch.obs=mvrnorm(batch.n, rep(0,p), (batch_sd_obs*diag(sds.lxtrue.train))^2)
      batch.temp=exp(batch_sort_train*batch.obs)
      # partially sort batch.obs for each gene
      batch.obs.sorted=c()
      for(j in 1:p){
        batch.id.unique2=batch.id.unique
        batch.id.unique.new=rep(NA, batch.n)
        batch.tempj=batch.temp[,j]
        for(b in 1:(batch.n-1)){
          probs=batch.tempj/sum(batch.tempj)
          sel=rmultinom(1, 1, probs)
          batch.id.unique.new[b]=batch.id.unique2[sel==1]
          batch.id.unique2=batch.id.unique2[sel!=1]
          batch.tempj=batch.tempj[sel!=1]
        }
        batch.id.unique.new[batch.n]=batch.id.unique[!(batch.id.unique %in% batch.id.unique.new)]
        batch.obsj.sorted=batch.obs[batch.id.unique.new,j]   
        batch.obs.sorted=cbind(batch.obs.sorted, batch.obsj.sorted)
      }
      batch.obs.mat=batch.obs.sorted[rep(1:batch.n, each=batch_size),]
      
      if(beta_sort_train==0){
        gamma.train.batch=gamma.train*exp(batch.mu.mat) ## gamma dist multiplied by a constant is still a gamma dist with scale par multiplied by the constant, which affects both mean and variance of the gamma dist
        id.new=id
      }else{
        temp=exp(beta_sort_train*log(t))
        id2=id
        id.new=rep(NA, n_train)
        for(i in 1:(n_train-1)){
          probs=temp/sum(temp)
          sel=rmultinom(1, 1, probs)
          id.new[i]=id2[sel==1]
          id2=id2[sel!=1]
          temp=temp[sel!=1]
        }
        id.new[n_train]=id[!(id %in% id.new)]
        gamma.train.batch=gamma.train[id.new,]*exp(batch.mu.mat)
      }
      gamma.train.batch.vec=c(gamma.train.batch)      
      xobs.vec=rpois(n_train*p, gamma.train.batch.vec*libsizes.sample.vec)
      batch.obs.vec=c(batch.obs.mat)
      xobs.vec=xobs.vec*exp(batch.obs.vec)  ## introduce batch effects on realized negative binomial counts
      x.train=matrix(log2(xobs.vec+0.5), nrow=n_train)  ## log2 scale
      row.sums.x.train=rowSums(x.train)
      allNA=1*(-p %in% row.sums.x.train)  ## indicator of existence of subjects whose all gene expressions are 0
      while(allNA==1){
        xobs.vec=rpois(n_train*p, gamma.train.batch.vec*libsizes.sample.vec)
        xobs.vec=xobs.vec*exp(batch.obs.vec)
        x.train=matrix(log2(xobs.vec+0.5), nrow=n_train)
        row.sums.x.train=rowSums(x.train)
        allNA=1*(-p %in% row.sums.x.train)
      }
      x.train.count=matrix(xobs.vec, nrow=n_train)      
    }else{
      id.new=id
      gamma.train.vec=c(gamma.train)
      xobs.vec=rpois(n_train*p, gamma.train.vec*libsizes.sample.vec)
      x.train=matrix(log2(xobs.vec+0.5), nrow=n_train)  ## log2 scale
      row.sums.x.train=rowSums(x.train)
      allNA=1*(-p %in% row.sums.x.train)  ## indicator of existence of subjects whose all gene expressions are 0
      while(allNA==1){
        xobs.vec=rpois(n_train*p, gamma.train.vec*libsizes.sample.vec)
        x.train=matrix(log2(xobs.vec+0.5), nrow=n_train)
        row.sums.x.train=rowSums(x.train)
        allNA=1*(-p %in% row.sums.x.train)
      }
      x.train.count=matrix(xobs.vec, nrow=n_train)   
    }

    
    ## test data

    # generate gamma-distributed true gene expression 
    norm.test=mvrnorm(n_test, mu=rep(0,p), Sigma=corr)
    unif.test=pnorm(norm.test)   
    gamma.test=matrix(NA, nrow=n_test, ncol=p)
    for(g in 1:p){
      gamma.test[,g]=qgamma(unif.test[,g], shape=1/disps.gene[g], scale=exp(betas.gene[g])*disps.gene[g])
    }
    
    # generate survival time based on gamma-distributed true gene expression
    xtrue.test=gamma.test*2*10^6  ## scale to observed counts
    lxtrue.test=log2(xtrue.test+0.5)
    sds.lxtrue.test=apply(lxtrue.test, 2, sd)
    h=h0*exp(lxtrue.test%*%b0)
    t0=rexp(n_test, h)
    ct=runif(n_test, 0, 300) 
    t.test=ifelse(t0<=ct, t0, ct)
    delta.test=1*(t0<=ct)  ## about 75% event rate
    
    id.test=1:n_test
    ncase.test=sum(delta.test)
    
    # sample library size (i.e. sequencing depth) from the observed library sizes of the original MSK sarcoma dataset
    temp.sample=sample(libsizes, n_test, replace=T)
    libsizes.sample=temp.sample+round(runif(n_test, -10000, 10000))  # add some perturbation to avoid identical library sizes
    libsizes.sample.vec=rep(libsizes.sample, p)  
    # libsizes.sample.vec=rep(10^6, n_test*p) 
    
    # introduce batch effect
    batch.id.test=rep(1:(n_test/batch_size+1), each=batch_size)[1:n_test]
    batch.id.test.unique=unique(batch.id.test)
    if(he_test==1){
      batch.n=length(batch.id.test.unique)
      batch.mu=mvrnorm(batch.n, rep(0,p), (batch_sd_mu*diag(sds.lxtrue.test))^2)  ## batch sd is proportional to gene sd
      batch.mu=apply(batch.mu, 2, sort, decreasing=T)  ## sort batch effect for each gene
      batch.mu.mat=batch.mu[rep(1:batch.n, each=batch_size),]
      
      batch.obs=mvrnorm(batch.n, rep(0,p), (batch_sd_obs*diag(sds.lxtrue.test))^2)
      batch.temp=exp(batch_sort_test*batch.obs)
      # partially sort batch.obs for each gene
      batch.obs.sorted=c()
      for(j in 1:p){
        batch.id.unique2=batch.id.test.unique
        batch.id.unique.new=rep(NA, batch.n)
        batch.tempj=batch.temp[,j]
        for(b in 1:(batch.n-1)){
          probs=batch.tempj/sum(batch.tempj)
          sel=rmultinom(1, 1, probs)
          batch.id.unique.new[b]=batch.id.unique2[sel==1]
          batch.id.unique2=batch.id.unique2[sel!=1]
          batch.tempj=batch.tempj[sel!=1]
        }
        batch.id.unique.new[batch.n]=batch.id.test.unique[!(batch.id.test.unique %in% batch.id.unique.new)]
        batch.obsj.sorted=batch.obs[batch.id.unique.new,j]   
        batch.obs.sorted=cbind(batch.obs.sorted, batch.obsj.sorted)
      }
      batch.obs.mat=batch.obs.sorted[rep(1:batch.n, each=batch_size),]
      
      if(beta_sort_test==0){
        gamma.test.batch=gamma.test*exp(batch.mu.mat)     
        id.new.test=id.test
      }else{
        temp=exp(beta_sort_test*log(t.test))
        id2=id.test
        id.new.test=rep(NA, n_test)
        for(i in 1:(n_test-1)){
          probs=temp/sum(temp)
          sel=rmultinom(1, 1, probs)
          id.new.test[i]=id2[sel==1]
          id2=id2[sel!=1]
          temp=temp[sel!=1]
        }
        id.new.test[n_test]=id.test[!(id.test %in% id.new.test)]
        gamma.test.batch=gamma.test[id.new.test,]*exp(batch.mu.mat)
      }
      gamma.test.batch.vec=c(gamma.test.batch)
      xobs.vec=rpois(n_test*p, gamma.test.batch.vec*libsizes.sample.vec)
      batch.obs.vec=c(batch.obs.mat)
      xobs.vec=xobs.vec*exp(batch.obs.vec) 
      x.test=matrix(log2(xobs.vec+0.5), nrow=n_test)  ## log2 scale
      row.sums.x.test=rowSums(x.test)
      allNA=1*(-p %in% row.sums.x.test)  ## indicator of existence of subjects whose all gene expressions are 0
      while(allNA==1){
        xobs.vec=rpois(n_test*p, gamma.test.batch.vec*libsizes.sample.vec)
        xobs.vec=xobs.vec*exp(batch.obs.vec)
        x.test=matrix(log2(xobs.vec+0.5), nrow=n_test)
        row.sums.x.test=rowSums(x.test)
        allNA=1*(-p %in% row.sums.x.test)
      }
      x.test.count=matrix(xobs.vec, nrow=n_test)      
    }else{
      id.new.test=id.test
      gamma.test.vec=c(gamma.test)
      xobs.vec=rpois(n_test*p, gamma.test.vec*libsizes.sample.vec)
      x.test=matrix(log2(xobs.vec+0.5), nrow=n_test)  ## log2 scale
      row.sums.x.test=rowSums(x.test)
      allNA=1*(-p %in% row.sums.x.test)  ## indicator of existence of subjects whose all gene expressions are 0
      while(allNA==1){
        xobs.vec=rpois(n_test*p, gamma.test.vec*libsizes.sample.vec)
        x.test=matrix(log2(xobs.vec+0.5), nrow=n_test)
        row.sums.x.test=rowSums(x.test)
        allNA=1*(-p %in% row.sums.x.test)
      }
      x.test.count=matrix(xobs.vec, nrow=n_test)   
    }

    
    ####### normalization ######

    if(norm_type==1){ # total count
      x.train_tc=norm.TC(raw=t(x.train.count))
      x.train.count=t(x.train_tc$dataNormalized)
      x.train=log2(x.train.count+0.5)

      #x.test_tc=norm.TCfrozen(raw=t(x.test),scalingFactor=x.train_tc$scalingFactor)
      x.test_tc=norm.TC(raw=t(x.test.count))
      x.test.count=t(x.test_tc$dataNormalized)
      x.test=log2(x.test.count+0.5)
      #x.test=t(t(x.test)/x.train_tc$scalingFactor)
    }

    if(norm_type==2){ # upper quantile
      x.train_uq=norm.UQ(raw=t(x.train.count))
      x.train.count=t(x.train_uq$dataNormalized)
      x.train=log2(x.train.count+0.5)

      x.test_uq=norm.TC(raw=t(x.test.count))
      x.test.count=t(x.test_uq$dataNormalized)
      x.test=log2(x.test.count+0.5)
    }

    if(norm_type==3){ # tmm
      x.train_tmm=norm.TMM(raw=t(x.train.count))
      x.train.count=t(x.train_tmm$dataNormalized)
      x.train=log2(x.train.count+0.5)

      x.test_tmm=norm.TMM(raw=t(x.test.count))
      x.test.count=t(x.test_tmm$dataNormalized)
      x.test=log2(x.test.count+0.5)
    }

    # if(norm_type==4){ # POISSONSEQ
    #   x.train_PoissonSeq=norm.PoissonSeq(raw=t(2^x.train))
    #   x.train=log(t(x.train_PoissonSeq$dataNormalized))
    #
    #   x.test_PoissonSeq=norm.PoissonSeq(raw=t(2^x.test))
    #   x.test=log(t(x.test_PoissonSeq$dataNormalized))
    # }

    if(norm_type==4){ # quantile normalization
      x.train=t(normalize.quantiles(t(x.train)))
      x.train.count=2^x.train
      x_quantile=sort(x.train[1,])

      x.test=t(apply(x.test, 1, function(u){x_quantile[rank(u)]}))
      x.test.count=2^x.test
    }

    colnames(x.train)=geneid
    colnames(x.test)=geneid


    #### pre-screen

    means_train=colMeans(x.train)
    x_train_sub=x.train[,means_train>0]
    xcorr_train=cor(x_train_sub)
    exclude_train=c()
    for(i in 1:ncol(x_train_sub)){
      for(j in i:ncol(x_train_sub)){
        if(xcorr_train[i,j]>ps_rho_cutoff & xcorr_train[i,j]!=1){
          exclude_train_1=ifelse(means_train[geneid==rownames(xcorr_train)[i]]>=means_train[geneid==colnames(xcorr_train)[j]],
                         colnames(xcorr_train)[j], rownames(xcorr_train)[i])
          exclude_train=c(exclude_train, exclude_train_1)
        }
      }
    }
    if(length(unique(exclude_train))==0){
      x_train=x_train_sub
    }else{
      x_train=x_train_sub[,!(colnames(x_train_sub) %in% unique(exclude_train))]
    }
    p_sub=ncol(x_train)
    subsize[iter]=p_sub

    x0_train=x.train[,nonzero_position]
    x0_test=x.test[,nonzero_position]

    #### gene standardization

    if(standardize==1){
      train.means=colMeans(x_train)
      train.sds=apply(x_train, 2, sd)
      x_train=(x_train-rep(1,nrow(x_train))%*%t(train.means))/(rep(1,nrow(x_train))%*%t(train.sds))

      train0.means=colMeans(x0_train)
      train0.sds=apply(x0_train, 2, sd)
      x0_train=(x0_train-rep(1,nrow(x0_train))%*%t(train0.means))/(rep(1,nrow(x0_train))%*%t(train0.sds))

      test.means=colMeans(x.test)
      test.sds=apply(x.test, 2, sd)
      x.test=(x.test-rep(1,nrow(x.test))%*%t(test.means))/(rep(1,nrow(x.test))%*%t(test.sds))

      test0.means=colMeans(x0_test)
      test0.sds=apply(x0_test, 2, sd)
      x0_test=(x0_test-rep(1,nrow(x0_test))%*%t(test0.means))/(rep(1,nrow(x0_test))%*%t(test0.sds))
    }


    data_train=as.data.frame(cbind(id.new, t[id.new], delta[id.new], batch.id, x_train))
    names(data_train)[1:3]=c("id","t","delta")
    
    data_test=as.data.frame(cbind(id.new.test, t.test[id.new.test], delta.test[id.new.test], batch.id.test, x.test))
    names(data_test)[1:4]=c("id","t","delta","batch.id")

    x_test=x.test
    

    ################ analysis ##############

    if(length(nonzero_position)>2){
      # oracle analysis
      coxfit_o=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x0_train)), error=function(e) e, warning=function(w) w)
      if(is(coxfit_o, "warning") | is(coxfit_o, "error")){
        coxfit_o=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x0_train), theta=1))
      }
      b_o0=rep(0,p)
      b_o0[nonzero_position]=summary(coxfit_o)$coefficient[,1]
      b_o[iter,]=b_o0
  
      coxfit_o_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test), init=summary(coxfit_o)$coefficient[,1], control=coxph.control(iter.max=0))
      c_o_test[iter]=summary(coxfit_o_test)$concordance[1]      
    }

    # univariate analysis with handling effect
    ps_u=lkhd_u=rep(0,p_sub)
    for(i in 1:p_sub){
      coxfit_ui=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,i]))
      if(!is.na(coef(coxfit_ui)) & abs(coef(coxfit_ui))<20){
        ps_u[i]=summary(coxfit_ui)$coefficient[1,5]
        lkhd_u[i]=coxfit_ui$loglik[2]
      }else{
        ps_u[i]=0.99
        lkhd_u[i]=-10000
      }
    }

    ps_u_sort=sort(ps_u)
    # cutoff_grid=seq(univ_cutoff_range[1], univ_cutoff_range[2], length.out=nuniv_cutoff)
    aics=rep(NA, round(min(ncase.train/10,p_sub/4)))
    for(k in 1:round(min(ncase.train/10,p_sub/4))){
      cut=ps_u_sort[k]
      selected=1*(ps_u<=cut)
      if(sum(selected, na.rm=T)!=0){
        coxfit_u=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,selected==1])),
                          error=function(e) e, warning=function(w) w)
        if(!(is(coxfit_u, "warning") | is(coxfit_u, "error"))){
          # if(max(abs(summary(coxfit_u)$coefficient[,1]))<15){
            aics[k]=extractAIC(coxfit_u)[2]
          # }
        }
      }
    }

    cut_sel=ps_u_sort[which(aics==min(aics, na.rm=T))[1]]
    sel_genes=colnames(x_train)[ps_u<=cut_sel]
    coxfit_u=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,ps_u<=cut_sel]))
    b_u0=rep(0,p)
    b_u0[geneid %in% sel_genes]=summary(coxfit_u)$coefficient[,1]
    b_u[iter,]=b_u0
    coxfit_u_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test[,geneid %in% sel_genes]), init=summary(coxfit_u)$coefficient[,1], control=coxph.control(iter.max=0))
    c_u_test[iter]=summary(coxfit_u_test)$concordance[1]


    lkhd_p_sort=sort(lkhd_u, decreasing=TRUE)
    lkhd_p_thres=lkhd_p_sort[round(min(ncase.train/10,p_sub/4))]
    sel_p_genes=colnames(x_train)[lkhd_u>=lkhd_p_thres]
    x_p_train=as.matrix(x_train[,lkhd_u>=lkhd_p_thres])
    inifit_p=tryCatch(coxph(Surv(data_train$t, data_train$delta)~as.matrix(x_train[,lkhd_u>=lkhd_p_thres])), error=function(e) e, warning=function(w) w)
    if(is(inifit_p, "warning") | is(inifit_p, "error")){
      inifit_p=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), theta=1))
    }else{
      if(max(abs(coef(inifit_p)))>10){
        inifit_p=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), theta=1))
      }
    }
    x_p_test=as.matrix(x_test[,geneid %in% sel_p_genes])

    print("mark1")

    # Lasso-penalized analysis with handling effect
    alpha_grid=seq(1,1,0.1)
    cvs_alpha=lambdas_alpha=c()
    for(alpha in alpha_grid){
      cv_l=cv.glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F)
      lambda_min=cv_l$lambda.min
      cvm_min=cv_l$cvm[cv_l$lambda==lambda_min]
      cvs_alpha=c(cvs_alpha, cvm_min)
      lambdas_alpha=c(lambdas_alpha, lambda_min)
    }
    lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    glmnet_l=glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F)
    b_l0temp=as.vector(glmnet_l$beta)
    b_l0=rep(0,p)
    b_l0[geneid %in% sel_p_genes]=b_l0temp
    b_l[iter,]=b_l0
    l_l[iter]=lambda_sel

    coxfit_l_test=coxph(Surv(data_test$t,data_test$delta)~x_p_test, init=b_l0temp, control=coxph.control(iter.max=0))
    c_l_test[iter]=summary(coxfit_l_test)$concordance[1]


    # Adaptive Lasso-penalized analysis with handling effect
    w_a=1/abs(coef(inifit_p))

    cvs_alpha=lambdas_alpha=c()
    for(alpha in alpha_grid){
      cv_a=cv.glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F, penalty.factor=w_a)
      lambda_min=cv_a$lambda.min
      cvm_min=cv_a$cvm[cv_a$lambda==lambda_min]
      cvs_alpha=c(cvs_alpha, cvm_min)
      lambdas_alpha=c(lambdas_alpha, lambda_min)
    }
    lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    glmnet_a=glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F, penalty.factor=w_a)
    b_a0temp=as.vector(glmnet_a$beta)
    b_a0=rep(0,p)
    b_a0[geneid %in% sel_p_genes]=b_a0temp
    b_a[iter,]=b_a0
    l_a[iter]=lambda_sel

    coxfit_a_test=coxph(Surv(data_test$t,data_test$delta)~x_p_test, init=b_a0temp, control=coxph.control(iter.max=0))
    c_a_test[iter]=summary(coxfit_a_test)$concordance[1]

    print("mark2")

    ################ stratified analysis ##############

    if(stratify_train==1){
      
      if(length(nonzero_position)>2){
      # oracle analysis
        coxfit_os=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x0_train)+strata(data_train$batch.id)),
                           error=function(e) e, warning=function(w) w)
        if(is(coxfit_os, "warning") | is(coxfit_os, "error")){
          coxfit_os=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x0_train), theta=1)+strata(data_train$batch.id))
        }
        b_os0=rep(0,p)
        b_os0[nonzero_position]=summary(coxfit_os)$coefficient[,1]
        b_os[iter,]=b_os0
  
        if(stratify_test==0){
          coxfit_os_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test), init=summary(coxfit_os)$coefficient[,1], control=coxph.control(iter.max=0))
        }
        if(stratify_test==1){
          coxfit_os_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test)+strata(data_test$batch.id), init=summary(coxfit_os)$coefficient[,1], control=coxph.control(iter.max=0))
        }
        c_os_test[iter]=summary(coxfit_os_test)$concordance[1]
        print("mark2-1")
      }
      
      # univariate analysis
      ps_us=lkhd_us=rep(0,p_sub)
      for(i in 1:p_sub){
        coxfit_usi=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,i]) + strata(data_train$batch.id)), error=function(e) e)
        if(!is(coxfit_usi, "error")){
          if(!is.na(coef(coxfit_usi))){
            if(abs(coef(coxfit_usi))<20){
              ps_us[i]=summary(coxfit_usi)$coefficient[1,5]
              lkhd_us[i]=coxfit_usi$loglik[2]
            }else{
              ps_us[i]=0.99
              lkhd_us[i]=-10000
            }
          }else{
            ps_us[i]=0.99
            lkhd_us[i]=-10000
          }
        }else{
          ps_us[i]=0.99
          lkhd_us[i]=-10000
        }
      }
      print("mark2-2")
      ps_us_sort=sort(ps_us)
      # cutoff_grid=seq(univ_cutoff_range[1], univ_cutoff_range[2], length.out=nuniv_cutoff)
      aics=rep(NA, round(min(ncase.train/10,p_sub/4)))
      for(k in 1:round(min(ncase.train/10,p_sub/4))){
        cut=ps_us_sort[k]
        selected=1*(ps_us<=cut)
        if(sum(selected, na.rm=T)!=0){
          coxfit_us=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,selected==1 & !is.na(selected)])+strata(data_train$batch.id)),
                            error=function(e) e, warning=function(w) w)
          if(!(is(coxfit_us, "warning") | is(coxfit_us, "error"))){
            # if(max(abs(summary(coxfit_us)$coefficient[,1]))<15){
              aics[k]=extractAIC(coxfit_us)[2]
            # }
          }
        }
      }
      print("mark2-3")
      cut_sel=ps_us_sort[which(aics==min(aics, na.rm=T))[1]]
      sel_genes_strat=colnames(x_train)[ps_us<=cut_sel]
      coxfit_us=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,ps_us<=cut_sel])+strata(data_train$batch.id))
      b_us0=rep(0,p)
      b_us0[geneid %in% sel_genes_strat]=summary(coxfit_us)$coefficient[,1]
      b_us[iter,]=b_us0

      if(stratify_test==0){
        coxfit_us_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test[,geneid %in% sel_genes_strat]), init=summary(coxfit_us)$coefficient[,1], control=coxph.control(iter.max=0))
      }
      if(stratify_test==1){
        coxfit_us_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test[,geneid %in% sel_genes_strat]) + strata(data_test$batch.id), init=summary(coxfit_us)$coefficient[,1], control=coxph.control(iter.max=0))
      }
      c_us_test[iter]=summary(coxfit_us_test)$concordance[1]
      print("mark2-4")

      lkhd_ps_sort=sort(lkhd_us, decreasing=TRUE)
      lkhd_ps_thres=lkhd_ps_sort[round(min(ncase.train/10,p_sub/4))]
      sel_p_genes_strat=colnames(x_train)[lkhd_us>=lkhd_ps_thres]
      x_ps_train=as.matrix(x_train[,lkhd_us>=lkhd_ps_thres])
      inifit_ps=tryCatch(coxph(Surv(data_train$t, data_train$delta)~as.matrix(x_train[,lkhd_us>=lkhd_ps_thres])+strata(data_train$batch.id)),
                         error=function(e) e, warning=function(w) w)
      if(is(inifit_ps, "warning") | is(inifit_ps, "error")){
        inifit_ps=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_us>=lkhd_ps_thres]), theta=1)+strata(data_train$batch.id))
      }else{
        if(max(abs(coef(inifit_ps)))>10){
          inifit_ps=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_us>=lkhd_ps_thres]), theta=1)+strata(data_train$batch.id))
        }
      }

      x_nps_test=as.matrix(x_test[,geneid %in% sel_p_genes_strat])

      print("mark3")

      # Lasso-penalized analysis with handling effect

      pfit_ls=glmnet_strat_cv(data=cbind(data_train[,1:4],x_train[,lkhd_us>=lkhd_ps_thres]), lambda_max=lambda_max_glmnet, nlambda=nlambda, nfold=nfold, penalty_wt=rep(1,length(coef(inifit_ps))), threshold=10^(-5))
      b_ls0=rep(0,p)
      b_ls0temp=pfit_ls[[1]]
      b_ls0[geneid %in% sel_p_genes_strat]=b_ls0temp
      b_ls[iter,]=b_ls0
      l_ls[iter]=pfit_ls[[2]]
      coxfit_ls=coxph(Surv(data_train$t,data_train$delta)~x_ps_train, init=b_ls0temp, control=coxph.control(iter.max=0))

      if(stratify_test==0){
        coxfit_ls_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test, init=b_ls0temp, control=coxph.control(iter.max=0))
      }
      if(stratify_test==1){
        coxfit_ls_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test + strata(data_test$batch.id), init=b_ls0temp, control=coxph.control(iter.max=0))
      }
      c_ls_test[iter]=summary(coxfit_ls_test)$concordance[1]

      # Adaptive Lasso-penalized analysis with handling effect
      w_a=1/abs(coef(inifit_ps))

      pfit_as=glmnet_strat_cv(data=cbind(data_train[,1:4],x_train[,lkhd_us>=lkhd_ps_thres]), lambda_max=lambda_max_glmnet, nlambda=nlambda, nfold=nfold, penalty_wt=w_a, threshold=10^(-5))
      b_as0=rep(0,p)
      b_as0temp=pfit_as[[1]]
      b_as0[geneid %in% sel_p_genes_strat]=b_as0temp
      b_as[iter,]=b_as0
      l_as[iter]=pfit_as[[2]]
      al_as[iter]=pfit_as[[3]]
      coxfit_as=coxph(Surv(data_train$t,data_train$delta)~x_ps_train, init=b_as0temp, control=coxph.control(iter.max=0))

      if(stratify_test==0){
        coxfit_as_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test, init=b_as0temp, control=coxph.control(iter.max=0))
      }
      if(stratify_test==1){
        coxfit_as_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test + strata(data_test$batch.id), init=b_as0temp, control=coxph.control(iter.max=0))
      }
      c_as_test[iter]=summary(coxfit_as_test)$concordance[1]
    }

    print("mark4")

    ################ analysis with ComBat ##############

    ####### ComBat ######

    # if(cmbt_train==1){
    #   if(cmbt_frozen==0){
    #     x_train_full=t(ComBat(dat=t(x.train.count), batch=data_train$batch.id, mod=NULL, par.prior=F, mean.only=T))
    #   }
    #   if(cmbt_frozen==1){
    #     temp=t(ComBat_v2(dat=t(x.train.count), batch=data_train$batch.id, mod=NULL, par.prior=F, mean.only=T, output.grand.mean.var.pooled=T))
    #     x_train_full=t(temp[[3]])
    #     grand.mean=temp[[1]]
    #     var.pooled=temp[[2]]
    #   }
    #   x_train_full_col_missing=colSums(is.na(x_train_full))
    #
    #   if(cmbt_test==1){
    #     if(cmbt_frozen==0){
    #       x_test=t(ComBat(dat=t(x.test.count), batch=data_test$batch.id, mod=NULL, par.prior=F, mean.only=T))
    #       x_test_col_missing=colSums(is.na(x_test))
    #     }
    #     if(cmbt_frozen==1){
    #       x_test=t(ComBat_v2(dat=t(x.test.count), batch=data_test$batch.id, mod=NULL, par.prior=F, mean.only=T, grand.mean=grand.mean, var.pooled=var.pooled))
    #       x_test_col_missing=colSums(is.na(x_test))
    #     }
    #   }else{
    #     x_test_col_missing=rep(0, ncol(x.test.count))
    #     x_test=x.test.count
    #   }
    #
    #   x.train=log(x_train_full[,x_train_full_col_missing==0 & x_test_col_missing==0])
    #   x.test=log(x_test[,x_train_full_col_missing==0 & x_test_col_missing==0])
    #   geneid2=geneid[x_train_full_col_missing==0 & x_test_col_missing==0]
    #   colnames(x.train)=geneid2
    #   colnames(x.test)=geneid2


    if(cmbt_train==1){
      if(cmbt_frozen==0){
        x_train_full=t(ComBat_seq(counts=t(x.train.count), batch=data_train$batch.id))
      }
      if(cmbt_frozen==1){
        temp=t(ComBat_v2(dat=t(x.train.count), batch=data_train$batch.id, mod=NULL, par.prior=F, mean.only=T, output.grand.mean.var.pooled=T))
        x_train_full=t(temp[[3]])
        grand.mean=temp[[1]]
        var.pooled=temp[[2]]
      }

      print("mark5")
      
      if(cmbt_test==1){
        if(cmbt_frozen==0){
          x_test=t(ComBat_seq(counts=t(x.test.count), batch=data_test$batch.id))
        }
        if(cmbt_frozen==1){
          x_test=t(ComBat_v2(dat=t(x.test.count), batch=data_test$batch.id, mod=NULL, par.prior=F, mean.only=T, grand.mean=grand.mean, var.pooled=var.pooled))
        }
      }else{
        x_test=x.test.count
      }

      print("mark6")
      
      x.train=log2(x_train_full+0.5)
      x.test=log2(x_test+0.5)
      colnames(x.train)=geneid
      colnames(x.test)=geneid

      #### pre-screen

      means_train=colMeans(x.train)
      x_train_sub=x.train[,means_train>0]
      xcorr_train=cor(x_train_sub)
      exclude_train=c()
      for(i in 1:ncol(x_train_sub)){
        for(j in i:ncol(x_train_sub)){
          if(xcorr_train[i,j]>ps_rho_cutoff & xcorr_train[i,j]!=1){
            exclude_train_1=ifelse(means_train[geneid==rownames(xcorr_train)[i]]>=means_train[geneid==colnames(xcorr_train)[j]],
                                   colnames(xcorr_train)[j], rownames(xcorr_train)[i])
            exclude_train=c(exclude_train, exclude_train_1)
          }
        }
      }
      if(length(unique(exclude_train))==0){
        x_train=x_train_sub
      }else{
        x_train=x_train_sub[,!(colnames(x_train_sub) %in% unique(exclude_train))]
      }
      p_sub=ncol(x_train)
      subsize_c[iter]=p_sub

      print("mark7")
      
      x0_train=x.train[,nonzero_position]
      x0_test=x.test[,nonzero_position]

      #### gene standardization

      if(standardize==1){
        train.means=colMeans(x_train)
        train.sds=apply(x_train, 2, sd)
        x_train=(x_train-rep(1,nrow(x_train))%*%t(train.means))/(rep(1,nrow(x_train))%*%t(train.sds))

        train0.means=colMeans(x0_train)
        train0.sds=apply(x0_train, 2, sd)
        x0_train=(x0_train-rep(1,nrow(x0_train))%*%t(train0.means))/(rep(1,nrow(x0_train))%*%t(train0.sds))

        test.means=colMeans(x.test)
        test.sds=apply(x.test, 2, sd)
        x.test=(x.test-rep(1,nrow(x.test))%*%t(test.means))/(rep(1,nrow(x.test))%*%t(test.sds))

        test0.means=colMeans(x0_test)
        test0.sds=apply(x0_test, 2, sd)
        x0_test=(x0_test-rep(1,nrow(x0_test))%*%t(test0.means))/(rep(1,nrow(x0_test))%*%t(test0.sds))
      }

      data_train=as.data.frame(cbind(id.new, t[id.new], delta[id.new], batch.id, x_train))
      names(data_train)[1:3]=c("id","t","delta")
      
      data_test=as.data.frame(cbind(id.new.test, t.test[id.new.test], delta.test[id.new.test], batch.id.test, x.test))
      names(data_test)[1:4]=c("id","t","delta","batch.id")
      
      x_test=x.test


      # oracle analysis

      if(length(nonzero_position)>2){
        coxfit_oc=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x0_train)),error=function(e) e, warning=function(w) w)
        if(is(coxfit_oc, "warning") | is(coxfit_oc, "error")){
          coxfit_oc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x0_train), theta=1))
        }
        b_oc0=rep(0,p)
        b_oc0[nonzero_position]=summary(coxfit_oc)$coefficient[,1]
        b_oc[iter,]=b_oc0
  
        coxfit_oc_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test), init=summary(coxfit_oc)$coefficient[,1], control=coxph.control(iter.max=0))
        c_oc_test[iter]=summary(coxfit_oc_test)$concordance[1]
      }

      # univariate analysis with handling effect
      ps_uc=lkhd_uc=rep(0,p_sub)
      for(i in 1:p_sub){
        coxfit_uci=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,i]))
        if(!is.na(coef(coxfit_uci)) & abs(coef(coxfit_uci))<20){
          ps_uc[i]=summary(coxfit_uci)$coefficient[1,5]
          lkhd_uc[i]=coxfit_uci$loglik[2]
        }else{
          ps_uc[i]=0.99
          lkhd_uc[i]=-10000
        }
      }

      ps_uc_sort=sort(ps_uc)
      # cutoff_grid=seq(univ_cutoff_range[1], univ_cutoff_range[2], length.out=nuniv_cutoff)
      aics=rep(NA, round(min(ncase.train/10,p_sub/4)))
      for(k in 1:round(min(ncase.train/10,p_sub/4))){
        cut=ps_uc_sort[k]
        selected=1*(ps_uc<=cut)
        if(sum(selected, na.rm=T)!=0){
          coxfit_uc=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,selected==1])),
                             error=function(e) e, warning=function(w) w)
          if(is(coxfit_uc, "warning") | is(coxfit_uc, "error")){
            coxfit_uc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,selected==1]), theta=1))
          }else{
            if(max(abs(coef(coxfit_uc)))>10){
              coxfit_uc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,selected==1]), theta=1))
            }
          }
          aics[k]=extractAIC(coxfit_uc)[2]
        }
      }

      cut_sel=ps_uc_sort[which(aics==min(aics, na.rm=T))[1]]
      sel_genes=colnames(x_train)[ps_uc<=cut_sel]
      coxfit_uc=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,ps_uc<=cut_sel]))
      if(max(abs(coef(coxfit_uc)))>10){
        coxfit_uc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,ps_uc<=cut_sel]), theta=1))
      }
      b_uc0=rep(0,p)
      b_uc0[geneid %in% sel_genes]=summary(coxfit_uc)$coefficient[,1]
      b_uc[iter,]=b_uc0

      coxfit_uc_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test[,geneid %in% sel_genes]), init=summary(coxfit_uc)$coefficient[,1], control=coxph.control(iter.max=0))
      c_uc_test[iter]=summary(coxfit_uc_test)$concordance[1]


      lkhd_pc_sort=sort(lkhd_uc, decreasing=TRUE)
      lkhd_pc_thres=lkhd_pc_sort[round(min(ncase.train/10,p_sub/4))]
      sel_pc_genes=colnames(x_train)[lkhd_uc>=lkhd_pc_thres]
      x_pc_train=as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres])
      inifit_pc=tryCatch(coxph(Surv(data_train$t, data_train$delta)~as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres])), error=function(e) e, warning=function(w) w)
      if(is(inifit_pc, "warning") | is(inifit_pc, "error")){
        inifit_pc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), theta=1))
      }else{
        if(max(abs(coef(inifit_pc)))>10){
          inifit_pc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), theta=1))
        }
      }
      x_pc_test=as.matrix(x_test[,geneid %in% sel_pc_genes])

      print("mark8")

      # Lasso-penalized analysis with handling effect
      alpha_grid=seq(1,1,0.1)
      cvs_alpha=lambdas_alpha=c()
      for(alpha in alpha_grid){
        cv_l=cv.glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F)
        lambda_min=cv_l$lambda.min
        cvm_min=cv_l$cvm[cv_l$lambda==lambda_min]
        cvs_alpha=c(cvs_alpha, cvm_min)
        lambdas_alpha=c(lambdas_alpha, lambda_min)
      }
      lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
      alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
      glmnet_lc=glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F)
      b_lc0temp=as.vector(glmnet_lc$beta)
      b_lc0=rep(0,p)
      b_lc0[geneid %in% sel_pc_genes]=b_lc0temp
      b_lc[iter,]=b_lc0
      l_lc[iter]=lambda_sel
      al_lc[iter]=alpha_sel
      coxfit_lc=coxph(Surv(data_train$t,data_train$delta)~x_pc_train, init=b_lc0temp, control=coxph.control(iter.max=0))
      coxfit_lc_test=coxph(Surv(data_test$t,data_test$delta)~x_pc_test, init=b_lc0temp, control=coxph.control(iter.max=0))
      c_lc_test[iter]=summary(coxfit_lc_test)$concordance[1]

      # Adaptive Lasso-penalized analysis with handling effect
      w_ac=1/abs(coef(inifit_pc))

      cvs_alpha=lambdas_alpha=c()
      for(alpha in alpha_grid){
        cv_a=cv.glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F, penalty.factor=w_ac)
        lambda_min=cv_a$lambda.min
        cvm_min=cv_a$cvm[cv_a$lambda==lambda_min]
        cvs_alpha=c(cvs_alpha, cvm_min)
        lambdas_alpha=c(lambdas_alpha, lambda_min)
      }
      lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
      alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
      glmnet_ac=glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F, penalty.factor=w_ac)
      b_ac0temp=as.vector(glmnet_ac$beta)
      b_ac0=rep(0,p)
      b_ac0[geneid %in% sel_pc_genes]=b_ac0temp
      b_ac[iter,]=b_ac0
      l_ac[iter]=lambda_sel
      al_ac[iter]=alpha_sel
      coxfit_ac=coxph(Surv(data_train$t,data_train$delta)~x_pc_train, init=b_ac0temp, control=coxph.control(iter.max=0))
      coxfit_ac_test=coxph(Surv(data_test$t,data_test$delta)~x_pc_test, init=b_ac0temp, control=coxph.control(iter.max=0))
      c_ac_test[iter]=summary(coxfit_ac_test)$concordance[1]
    }
    
    cat(paste("Iteration-",iter," ",sep=""))
    }, timeout=3600)},
    TimeoutException = function(ex) {message("Timeout. Skipping.")
      iter=iter-1
      print(paste("iter=",iter,sep=""))},
    # error=function(e){iter=iter-1}, warning=function(w){iter=iter-1})
    error=function(e){message("Error. Skipping.")
      iter=iter-1
      print(paste("iter=",iter,sep=""))}
    )  ## if there is any error or time out in the entire iteration, the iteration is repeated with a different seed until no error
  }
  cat(" done \n")
  
  
  for(second in c("o","u","l","a")){
    bname=paste("b_",second,sep="")
    bobj=eval(parse(text=bname))
    write.csv(bobj, paste(bname,"_",addon,"_",seed,".csv",sep=""), row.names=F)
    
    cname_test=paste("c_",second,"_test",sep="")
    cobj_test=eval(parse(text=cname_test))
    write.csv(cobj_test, paste(cname_test,"_",addon,"_",seed,".csv",sep=""), row.names=F)
    
    if(second %in% c("l","a")){
      lname=paste("l_",second,sep="")
      lobj=eval(parse(text=lname))
      write.csv(lobj, paste(lname,"_",addon,"_",seed,".csv",sep=""), row.names=F)      
    }
  }
  
  if(stratify_train==1){
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"s",sep="")
      bobj=eval(parse(text=bname))
      write.csv(bobj, paste(bname,"_",addon,"_",seed,".csv",sep=""), row.names=F)
      
      cname_test=paste("c_",second,"s_test",sep="")
      cobj_test=eval(parse(text=cname_test))
      write.csv(cobj_test, paste(cname_test,"_",addon,"_",seed,".csv",sep=""), row.names=F)
      
      if(second %in% c("l","a")){ 
        lname=paste("l_",second,"s",sep="")
        lobj=eval(parse(text=lname))
        write.csv(lobj, paste(lname,"_",addon,"_",seed,".csv",sep=""), row.names=F)
      }
    }
  }
  
  write.csv(subsize,paste("pssize_",addon,"_",seed,".csv",sep=""), row.names=F)
  
  if(cmbt_train==1){
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"c",sep="")
      bobj=eval(parse(text=bname))
      write.csv(bobj, paste(bname,"_",addon,"_",seed,".csv",sep=""), row.names=F)
      
      cname_test=paste("c_",second,"c_test",sep="")
      cobj_test=eval(parse(text=cname_test))
      write.csv(cobj_test, paste(cname_test,"_",addon,"_",seed,".csv",sep=""), row.names=F)
      
      if(second %in% c("l","a")){
        lname=paste("l_",second,"c",sep="")
        lobj=eval(parse(text=lname))
        write.csv(lobj, paste(lname,"_",addon,"_",seed,".csv",sep=""), row.names=F)
      }
    }
    write.csv(subsize_c,paste("pssizeC_",addon,"_",seed,".csv",sep=""), row.names=F)
  }
  }






