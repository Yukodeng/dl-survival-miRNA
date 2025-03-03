#### Deep Augmented miRNA-seq Data: applications of baseline survival models ####
# Date      Notes  
#   Feb2025 Issues: 1) augmented gene counts have negative values; gene number (p) does not align
#                   2) algorithm for introducing correlation between batch and survival (beta_sort)
#                   3) flag=1lambda=NA in penalized Lasso; "Error. Skipping" messages
# 03Mar2025 Changed log2 transformation to natural log

mainSimSeqedgeRpar.aug = function(
    beta0, nonzero_position, n_train, n_test, beta_sort_train=0, beta_sort_test=0, 
    train_batch_size=30, test_batch_size=20, n_batch=10, #maf augmented batch data
    he_train, he_test, batch_sd_mu, batch_sd_obs, batch_sort_train, batch_sort_test, stratify_train, stratify_test, cmbt_train, cmbt_test, cmbt_frozen,
    norm_type, norm_train, norm_test, standardize=0, ps_rho_cutoff=0.9, univ_cutoff_range=c(0,0.05), nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5, 
    addon, nsim, seed) {
  
  library(survival)
  library(preprocessCore)
  library(glmnet)
  library(sva)
  library(MASS)
  library(stringr)
  library(edgeR)
  library(R.utils)

  cat(paste0(addon,"-seed-",seed))
  
  options(stringsAsFactors=F)
  setwd("~/dl-survival-miRNA/simulation/batch_effect/")
  source("norm.R")
  source("glmnet_strat_cv_2.R")
  
  geneid=read.csv(file.path('augment-parametric simulation', "gene.names.538.csv"), header=T)$gene
  p=length(geneid)
  
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
  
  for(iter in (nsim*(seed-1)+1):(nsim*(seed-1)+nsim)){
    # print(iter)
    tryCatch(expr = {withTimeout({
      
      ##### Load data ###########
      train.dat=test.dat=data.frame(matrix(NA, nrow=0, ncol=p*2))
      for(i in 1:n_batch){
        ## train data
        train.sub=read.csv(
          file.path("~/dl-survival-miRNA/data/with_without_batch_maf", paste0("batch_",i,"_epochES_batch01_maf_generated",".csv")),
          header=T, skip=50*(iter-1), nrow=30)[,1:(p*2)]  ## train data n=30 from each batch
        
        ## test data
        test.sub=read.csv(
          file.path("~/dl-survival-miRNA/data/with_without_batch_maf", paste0("batch_",i,"_epochES_batch01_maf_generated",".csv")),
          header=T, skip=50*(iter-1)+30, nrow=20)[,1:(p*2)] ## test data n=20 from each batch
        
        ## merge batches
        colnames(train.sub)=colnames(test.sub)=rep(geneid,2)
        train.dat = rbind(train.dat, train.sub)
        test.dat = rbind(test.dat, test.sub)
      }
      
      ##### Note: keep the below 2 lines for now to prevent negative gene counts
      train.dat=abs(train.dat) 
      test.dat=abs(test.dat)
      
      # separate clean data and batch-contaminated data for train and test sets
      xtrue.train = train.dat[, 1:p]
      xtrue.test = test.dat[, 1:p]
      batch.train = train.dat[, (p+1):(p*2)]
      batch.test = test.dat[, (p+1):(p*2)]
      
      # get batch id
      batch.id=rep(1:n_batch, each=train_batch_size)[1:n_train]
      batch.id.test=rep(1:n_batch, each=test_batch_size)[1:n_test]
      
      
      # generate survival time based on gamma-distributed true gene expression
      ## train
      # 3/3/2025 (YD) Changed to natural logs per Andy
      # lxtrue.train=log2(xtrue.train+0.5)
      lxtrue.train=log(xtrue.train+0.5)
      sds.lxtrue.train=apply(lxtrue.train, 2, sd)
      # h0=1.5  ## for small gene effect size
      h0=2200  ## for large gene effect size
      b0=rep(0,p)
      b0[nonzero_position]=beta0
      h=h0*exp(as.matrix(lxtrue.train)%*%as.matrix(b0))
      t0=rexp(n_train, h)
      ct=runif(n_train, 0, 300) 
      t=ifelse(t0<=ct, t0, ct)
      delta=1*(t0<=ct)  ## about 75% event rate
      # sum(delta)/length(delta)
      id=id.new=1:n_train
      ncase.train=sum(delta)
      
      ## test
      # lxtrue.test=log2(xtrue.test+0.5)
      lxtrue.test=log(xtrue.test+0.5)
      sds.lxtrue.test=apply(lxtrue.test, 2, sd)
      h=h0*exp(as.matrix(lxtrue.test)%*%b0)
      t0=rexp(n_test, h)
      ct=runif(n_test, 0, 300) 
      t.test=ifelse(t0<=ct, t0, ct)
      delta.test=1*(t0<=ct)  ## about 75% event rate
      # sum(delta.test)/length(delta.test)
      id.test=id.new.test=1:n_test
      ncase.test=sum(delta.test)
      
      if(he_train==0) {x.train.count=xtrue.train} else {x.train.count=batch.train}
      if(he_test==0) {x.test.count=xtrue.test} else {x.test.count=batch.test}
      
      
      ##### Normalization #########
      
      if(norm_type==1){ # total count
        x.train_tc=norm.TC(raw=t(x.train.count))
        x.train.count=t(x.train_tc$dataNormalized)
        # x.train=log2(x.train.count+0.5)
        x.train=log(x.train.count+0.5)
        
        #x.test_tc=norm.TCfrozen(raw=t(x.test),scalingFactor=x.train_tc$scalingFactor)
        x.test_tc=norm.TC(raw=t(x.test.count))
        x.test.count=t(x.test_tc$dataNormalized)
        # x.test=log2(x.test.count+0.5)
        x.test=log(x.test.count+0.5)
        #x.test=t(t(x.test)/x.train_tc$scalingFactor)
      }
      
      if(norm_type==2){ # upper quantile
        x.train_uq=norm.UQ(raw=t(x.train.count))
        x.train.count=t(x.train_uq$dataNormalized)
        # x.train=log2(x.train.count+0.5)
        x.train=log(x.train.count+0.5)
        
        x.test_uq=norm.TC(raw=t(x.test.count))
        x.test.count=t(x.test_uq$dataNormalized)
        # x.test=log2(x.test.count+0.5)
        x.test=log(x.test.count+0.5)
      }
      
      if(norm_type==3){ # tmm
        x.train_tmm=norm.TMM(raw=t(x.train.count))
        x.train.count=t(x.train_tmm$dataNormalized)
        # x.train=log2(x.train.count+0.5)
        x.train=log(x.train.count+0.5)

        x.test_tmm=norm.TMM(raw=t(x.test.count))
        x.test.count=t(x.test_tmm$dataNormalized)
        # x.test=log2(x.test.count+0.5)
        x.test=log(x.test.count+0.5)
      }
      
      if(norm_type==4){ # quantile normalization
        x.train=t(normalize.quantiles(t(x.train)))
        # x.train.count=2^x.train
        x.train.count=exp(x.train)
        x_quantile=sort(x.train[1,])
        
        x.test=t(apply(x.test, 1, function(u){x_quantile[rank(u)]}))
        # x.test.count=2^x.test
        x.test.count=exp(x.test)
      }
      # colnames(x.train)=geneid
      # colnames(x.test)=geneid
      
      # pre-screen #
      
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
      
      # gene standardization #
      
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
      # View(data_test)
      x_test=x.test
      
      
      ##### Analysis ##############
      
      if(length(nonzero_position)>2){
        
        # oracle analysis -----------
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
      
      # univariate analysis with handling effect -----
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
      
      # Lasso-penalized analysis with handling effect ----
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
      
      
      # Adaptive Lasso-penalized analysis with handling effect -----
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
      
      
      # stratified analysis -----------
      
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
      
      # analysis with ComBat ----------------
      
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
        
        # x.train=log2(x_train_full+0.5)
        # x.test=log2(x_test+0.5)
        x.train=log(x_train_full+0.5)
        x.test=log(x_test+0.5) 
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
      cat(paste0("Iteration-",iter," "))
    }, timeout=7200) },
    TimeoutException = function(ex) { message("Timeout. Skipping.")
      iter=iter-1
      print(paste0("iter=",iter=))
    },
    # error=function(e){iter=iter-1}, warning=function(w){iter=iter-1})
    error=function(e){ message("Error. Skipping.")
      iter=iter-1
      print(paste0("iter=",iter))
    }
    ) ## if there is any error or time out in the entire iteration, the iteration is repeated with a different seed until no error
  }
  cat(" done \n")
  
  ##### save data ###########
  dir.create(file.path("results", addon), showWarnings=F)
  setwd(file.path("~/dl-survival-miRNA/simulation/batch_effect/results", addon))
  
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


# beta0=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)
# nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)
# n_train=300
# n_test=200
# ps_rho_cutoff=0.9
# beta_sort_train=0
# beta_sort_test=0
# train_batch_size=30
# test_batch_size=20
# n_batch=10 #maf augmented batch data
# he_train=0
# he_test=0
# batch_sd_mu=0
# batch_sd_obs=0
# batch_sort_train=5
# batch_sort_test=5
# stratify_train=1
# stratify_test=1
# cmbt_train=0
# cmbt_test=0
# cmbt_frozen=0
# norm_type=1
# norm_train=1
# norm_test=1
# standardize=0
# univ_cutoff_range=c(0,0.05)
# nuniv_cutoff=20
# lambda_max_glmnet=0.3
# nlambda=30
# nfold=5
# addon='he0000bs20norm00tsort00-aug'
# nsim=20
# seed=1
# 
# ##################### Large gene effects, mean-level batch effects only ###################################
# ### no train HE, no test HE, no norm, batch size 20, largest batch effect ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,0,0,0,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm00tsort00-aug',20,seed)
# 
# ### no train HE, no test HE, TC norm, batch size 20, largest batch effect ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,1,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm11tsort00-aug',20,seed)
# 
# # no train HE, no test HE, upper quantile norm, batch size 10 ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,2,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm22tsort00-aug',20,seed)
# 
# ### no train HE, no test HE, TMM norm, batch size 10 ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,3,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm33tsort00-aug',20,seed)
# 
# ### no train HE, no test HE, quantile norm, batch size 10 ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,4,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm44tsort00-aug',20,seed)
# 
# 
# ############################# train HE, no test HE #############################
# ### train HE, no test HE, no norm ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,0,0,0,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm00tsort00-aug',20,seed)
# 
# ### train HE, no test HE, TC norm, batch size 20 ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,1,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm11tsort00-aug',20,seed)
# 
# ### train HE, no test HE, upper quantile norm, batch size 20 ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,2,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm22tsort00-aug',20,seed)
# 
# ### train HE, no test HE, TMM norm, batch size 20 ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,3,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm33tsort00-aug',20,seed)
# 
# ### train HE, no test HE, quantile norm, batch size 20 ###
# mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,4,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm44tsort00-aug',20,seed)
# 
