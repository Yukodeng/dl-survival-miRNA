############################################################################# ##
# Date      Note
# 25Apr2025 The 10 true markers now are included in the 30 weak miRNA markers
# 02May2025 Cap augmented marker values by real-world data GLOBAL (max + sd) + noise
# 03Jun2025 Updates: 1) Vary test set across 20 iterations. 2) Add new quadratic 
#            with intercept survival simulation. 3) Add true nonlinear Cox model
# 12Jun2025 Updated augmented data preprocessing to cap values at log2 scale:
#            Cap = global max + 0.1*sd + e; e~Norm(0, 0.01*sd)
# 16Jun2025 Corrected columns from which to extract clean vs. batch-contaminate
#            augmented data (for version sent on 15Apr2025)
# 19Jun2025 Updated reticulate and slurm script generation to use dl-surv python 3.10 virtual env
# 16Jul2025 Rescaled augmented batch effect sizes by 1/2;
#           Added stratified Oracle (linear & nonlinear); pending Lasso implementation (runs slow)
############################################################################# ##

library(ggsurvfit)
library(survival)
library(edgeR)
library(glmnet)
library(dplyr)
library(ggplot2)
library(preprocessCore)
# library(splitstackshape) # for train-test split on multiple outcome variables

# Import Python package scikit-learn for train-test split
library(reticulate)
reticulate::use_python("~/dl-env/bin/python3.8", required = T)
skl_ms = reticulate::import("sklearn.model_selection")

options(stringsAsFactors=F)
Sys.setenv("RCPP_PARALLEL_BACKEND" = "single")
options(preprocessCore.nthreads = 1)

setwd("~/dl-survival-miRNA")


# Load Raw Augmented Gene Counts ---------------------------------------------

date = format(Sys.Date(), "%m%d%y")
N = 10000
P = 538 
n_batch = 10

test.ids = read.csv(file.path('data', 'test_ids_10000.txt'), sep='\t')$id # test set row indices
g.names = read.csv("data/augmentation_miRNA_names_538.csv")$gene      # gene names
nonzero_position = c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)   # RWD true gene positions
nonzero.genes = g.names[nonzero_position]


## Augmented data with/without batch (attempt APR112025): N=10000 each batch
### NOTE: Avoid repeated loading inside functions

cat("Percentages of Inf values in each batch:\n")

surv.data=data.frame(matrix(NA, nrow=0, ncol=P*2))
for(i in 1:n_batch) {
  ## read in each batch
  dat.sub=read.csv(
    file.path("raw-data", "with_without_batch_maf", paste0("batch-",i,"_epochES_maf_generated.csv")),
    header=F
  )[1:N,1:(P*2)]  
  
  ## percentage of inf values in each batch
  cat(glue::glue("Batch {i}:"), sum(is.infinite(unlist(dat.sub))) / (N*P) * 100,"%\n")

  ## merge batches
  colnames(dat.sub)=g.names
  surv.data = rbind(surv.data, dat.sub)
}



# Main Simulation Functions --------------------------------------------------

simulate_T <- function(h0, h.log, n_train, max.censor.time) {
  
  # simulate uncensored survival time
  t0 = -log(runif(n_train)) * h0 / exp(h.log)
  
  # create censoring time
  ct=runif(n_train, 0, max.censor.time) ## capped by a max time
  t=ifelse(t0<=ct, t0, ct)
  delta=1*(t0<=ct)
  ncase.train=sum(delta)
  
  hist(t) ## visualize (censored) survival time
  cat(glue::glue("Event rate: {ncase.train / n_train}\n\n")) ## about 75% event rate
  # print(summary(t)) ## summarize simulated survival time
  
  return(list(time=t, status=delta))
}


N=10000
N.test=10000
rwd.max=22.516867
rwd.sd=4.039547
test.ids=test.ids
sim.dataType="nl-shiftquad" #"nl-quadratic"
h0=19
max.censor.time=200
nonzero.genes=nonzero.genes
p=10
p.weak=30
plot_km=F
save_surv_data=F
save_beta=F
save_gene_data=F
n_batch=10
he_train=1
he_test=1
beta_sort_train=0.01
beta_sort_test=0
norm_type=0
run_analysis=F
stratify=1
test_size=1000
subset_sizes = c(100,500,1000,2000,5000,10000)
runs_per_size = c(20,20,20,20,20,20)
splits_per_size = c(3,5,5,10,10,10)
generate_subFile=T
GPUtype='l4'
seed=1234
sim.survdata <- function(surv.data,
                         N=10000,
                         N.test=10000,
                         rwd.max=22.516867,
                         rwd.sd=4.039547,
                         test.ids=NULL,
                         sim.dataType=NULL, # E.g., "linear-moderate", "nl-quadratic"
                         h0=5, 
                         max.censor.time=200, 
                         nonzero.genes=NULL,
                         p=10, 
                         p.weak=30,
                         plot_km=T,
                         save_surv_data=F,
                         save_beta=T,
                         save_gene_data=T,
                         n_batch=10,
                         he_train=0,
                         he_test=0,
                         # batch_sort_train=5,
                         # batch_sort_test=5,
                         beta_sort_train=0,
                         beta_sort_test=0,
                         norm_type=0,
                         run_analysis=T,
                         stratify=1, 
                         test_size=1000, 
                         subset_sizes = c(100,500,1000,2000,5000,10000),
                         runs_per_size = c(20,20,20,20,20,20),
                         splits_per_size = c(3,5,5,10,10,10),
                         generate_subFile=T, 
                         GPUtype='l4',
                         seed=1234) {
  
  set.seed(seed)
  
  source("simulation/norm.R")
  source("simulation/glmnet_strat_cv_2.R")
  norm_map = c('None','TC','UQ','TMM','Quantile','DEseq','Med')
  
  ## 1. Preprocess Gene Expression Data -------------------------
  
  cat("\n\nStarting survival simulation..\n\n")
  
  # Read in gene names
  g.names = read.csv("data/augmentation_miRNA_names_538.csv")$gene ## gene names
  P = length(g.names)
 
  ### Load augmented data with/without batch -----
  # Check if there's any negative values
  surv.data[surv.data<0] = 0
  
  # Cap extremely large values
  ## [update 6/12/2025] Cap each marker (on the log2 scale) by (max + 0.1*sd) 
  ## estimated from the 108 real-world sarcoma samples; then add Gaussian noise 
  ## to the log2-transformed counts with mean=0 and standard dev = 0.01*sd
  cap_and_add_noise <- function(x, cap, sd) {
    new_cap = cap + 0.1*sd
    xs = x[is.infinite(x) | (x > new_cap)]
    xs_hat = pmin( xs, new_cap + rnorm(length(xs), 0, 0.01*sd) )
    x[is.infinite(x) | (x > new_cap)] = xs_hat
    return(x)
  }
  
  surv.data.cap = 2^as.data.frame(lapply(log2(surv.data+1), cap_and_add_noise, cap=rwd.max, sd=rwd.sd)) - 1
  
  # Separate clean and batch-contaminated data
  surv.data.batch = surv.data.cap[, 1:P]
  surv.data.true = surv.data.cap[, (P+1):(P*2)]
  colnames(surv.data.true)=colnames(surv.data.batch)=g.names # rename gene ids in batch-present data
  

  ## 2. Select True Predictors From Augmented Data ------------

  workingt = t(surv.data.true)  
  
  ## Filter genes by count per million = 10
  cpm=cpm(workingt)
  keep=(rowSums(cpm>2)>=10) # adjustable
  workingt=workingt[keep,]
  
  ## Select genes with higher TC and variance
  tc=rowSums(workingt)
  mid.tc=(tc>median(tc) & tc<quantile(tc, 0.75)) # select genes with total count in the upper quartile
  cv=apply(workingt, 1, sd)/rowMeans(workingt)
  # update 050725: lowered quantile -> 0.4 so that there are enough qualifying genes
  mid.cv=(cv > quantile(cv, 0.4)) # select genes with higher variance (SD/mean)
  candidate.geneid=g.names[mid.tc & mid.cv]
  
  if (is.null(nonzero.genes)) {
    ### Select p (default=10) true genes associated with survival
    if (p<=length(candidate.geneid)) {
      selected.geneid = sample(candidate.geneid, p, replace=F)
    } else {
      selected.geneid = candidate.geneid
    }
  } else {
    ### Alternatively import true gene names
    selected.geneid = nonzero.genes
  }

  
  ## 3. Simulate Survival Outcome -------------------------

  # log-transform selected true genes
  xtrue = workingt[selected.geneid, ]
  xtrue.log = log(xtrue+0.5) # log2(xtrue+0.5)
  n_train = ncol(xtrue)
  
  ### 10 linear moderate signals --------
  if (sim.dataType=='linear-moderate') {
    
    cat("Simulating from linear risk function with moderate signals\n")
    
    ## calculate beta's as inverse std of gene expression
    sd.genes = apply(xtrue.log, 1, sd)
    beta0 = rep(c(1,-1),p/2) * round(2/sd.genes, 1)
    
    ## calculate log risk
    h.log = t(xtrue.log) %*% beta0
    
  ### 30 linear weak signals ------------
  } else if (sim.dataType=='linear-weak') {
      
    cat("Simulating from linear risk function with weak signals\n")

    ## select (default=30) true genes associated with survival time
    if (p.weak <= length(candidate.geneid)) {
      selected.geneid = 
        c(selected.geneid,
          sample(candidate.geneid[!candidate.geneid %in% selected.geneid], p.weak-p, replace=F)
          )
    } else {
      selected.geneid = candidate.geneid
    }
    # print(selected.geneid)
    
    # Log-transform gene expression
    xtrue.weak = workingt[selected.geneid,]
    # xtrue.weak.log = log2(xtrue.weak+0.5)
    xtrue.weak.log = log(xtrue.weak+0.5)
    
    # Simulate beta values
    sd.genes.weak = apply(xtrue.weak.log, 1, sd)
    ## Note: shrink beta coefficients to control the total signal quantity
    beta.scale = p/p.weak
    beta0 = rep(c(1,-1),p/2) * round(2/sd.genes.weak*beta.scale, 1)
    
    h.log = t(xtrue.weak.log) %*% beta0
  
  ### 10 quadratic terms -----------------
  } else if (sim.dataType=='nl-quadratic') {
    
    cat("Simulating with quadratic transformation\n")
    
    ## Scale marker expression data
    xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
    xtrue.nl.log = as.matrix(xtrue.log.scaled)^2
    
    sd.genes = apply(xtrue.nl.log, 1, sd)
    beta0 = rep(c(1,-1),5) * round(2/sd.genes, 1)
    
    h.log = t(xtrue.nl.log) %*% beta0
    
  #### 10 quadratic with intercept --------
  } else if (sim.dataType == 'nl-shiftquad') {
    
    cat("Simulating with (X - median)^2 transformation\n")

    ## Scale marker expression data
    xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
    c_val = apply(xtrue.log.scaled, 1, median) # centering values (median)
    xtrue.nl.log = sweep(xtrue.log.scaled, 1, c_val, "-")^2
    
    sd.genes = apply(xtrue.nl.log, 1, sd)
    beta0 = rep(c(1, -1), length.out=p) * round(2/sd.genes, 1)
    
    h.log = t(xtrue.nl.log) %*% beta0
    
  ### 10 interaction terms ---------------
  } else if (sim.dataType=='nl-interaction') {
    
    cat("Simulating with interaction effects\n")
   
    interact.gene = sample(selected.geneid, 1)
    cat(glue::glue("Interaction term: {interact.gene}\n\n"))
    
    xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
    interact.gene.val = xtrue.log.scaled[interact.gene,]
    
    xtrue.nl.log = sweep(xtrue.log.scaled, 2, interact.gene.val, FUN="*") #as.matrix(xtrue.log.scaled) %*% diag(interact.gene.val)
    sd.genes = apply(xtrue.nl.log, 1, sd)
    beta0 = rep(c(1,-1),5) * round(2/sd.genes, 1) 
    
    h.log = t(xtrue.nl.log) %*% beta0
    
  ### 10 sine terms ----------------------
  } else if (sim.dataType=='nl-sine') {

    cat("Simulating with non-linear sine-transformed terms\n\n")
    
    xtrue.nl.log = t(sin(t(xtrue.log))) # sine transformation
    
    sd.genes = apply(xtrue.nl.log, 1, sd)
    beta0 = rep(c(1,-1),5) * round(2/sd.genes, 2)
    
    h.log = t(xtrue.nl.log) %*% beta0
  }
  
  cat('Selected true markers:\n')
  print(selected.geneid)
  # set.seed(seed)
  
  ## Simulated survival time and censoring outcome
  surv.out = simulate_T(h0, h.log, n_train, max.censor.time)
  t = surv.out$time
  delta = surv.out$status
  
  
  ## 4. Introduce Batch Effects ----------------------------

  ### Split into train and test (stratified on events and batch) -----
  batch.id=rep(1:n_batch, each=N)
  batch.event.temp =
    as.data.frame(cbind(delta, batch.id)) %>%
    dplyr::rename(status=V1) %>%
    tibble::rownames_to_column('ID') %>%
    mutate(ID = as.numeric(ID))
  
  if (is.null(test.ids)) {
    test.ids = splitstackshape::stratified(
        batch.event.temp, c('status','batch.id'), size=N.test/nrow(batch.event.temp)
      ) %>% pull(ID)
  }

  xtrue.train = surv.data.true[-test.ids,] 
  xtrue.test  = surv.data.true[test.ids, ] 
  batch.train = surv.data.batch[-test.ids, ]
  batch.test  = surv.data.batch[test.ids, ]
  # table(batch.id[test.ids]) #check that test set has balanced batches
  
  t.train = t[-test.ids] 
  t.test = t[test.ids]
  delta.train = delta[-test.ids]
  delta.test = delta[test.ids]

  # reset index
  N.train = nrow(xtrue.train)
  rownames(xtrue.train)=rownames(batch.train)=names(t.train)=names(delta.train)=1:N.train
  rownames(xtrue.test)=rownames(batch.test)=names(t.test)=names(delta.test)=1:N.test
  
  ### HE Train ------------

  id=id.new=1:N.train
  batch.id.train = batch.id[-test.ids]
  batch.id.unique = unique(batch.id.train)

  if(he_train==0) {
    x.train.count=xtrue.train
    x.train=log(x.train.count+0.5)
    
  } else {
    xtrue.train.log = log(xtrue.train+0.5)
    batch.train.log = log(batch.train+0.5)
    ## NOTE: extract batch on the natural log scale
    ### update 7/16: scale batch effects by 1/2
    batch.obs = (batch.train.log - xtrue.train.log) / 2
    
    ### Note: skip partial sorting 
    batch.train.log = batch.obs + xtrue.train.log
    batch.train = exp(as.data.frame(lapply(batch.train.log, cap_and_add_noise, cap=rwd.max, sd=rwd.sd)))
    
    ### [Per Andy 7/21/2025] ----
    ### The batch effects should be sorted for all genes
    ### at the same time to preserve the correlation structure of the batch effects across genes. 
    ### You can calculate the median of the median batch effects across genes in each batch 
    ### and then simply “hard” sort the medians across batches. There is no need for partial sorting at this stage.
    
    ## get per batch median
    batch.obs$batch_id = batch.id.train
    batch.order = cbind(
      Median = apply(batch.obs %>% group_by(batch_id) %>% summarise(across(everything(), median)) %>% select(-batch_id), 1, median),
      batch_id = seq(1,10)
    ) %>% as.data.frame() %>% 
      arrange(desc(Median)) %>% pull(batch_id)
    
    batch.obs$batch_id = factor(batch.obs$batch_id, levels = batch.order)
    batch.obs.sorted = batch.obs[order(batch.obs$batch_id), ] %>% select(-batch_id)
    
    # batch.obs.med = batch.obs %>%
    #   group_by(batch_id) %>%
    #   summarise(across(everything(), median)) %>%
    #   select(-batch_id) %>%
    #   as.matrix()
    # ## partially sort batch.obs for each gene
    # batch.temp = exp(batch_sort_train*batch.obs.med)
    # batch.obs.sorted = matrix(nrow=N.train, ncol=P)
    # colnames(batch.obs.sorted) = g.names
    # 
    # for(j in 1:P){
    #   batch.id.unique2=batch.id.unique
    #   batch.id.unique.new=rep(NA, n_batch)
    #   batch.tempj=batch.temp[,j]
    #   
    #   for(b in 1:(n_batch-1)){
    #     probs = batch.tempj/sum(batch.tempj)
    #     sel = rmultinom(1, 1, probs)
    #     batch.id.unique.new[b] = batch.id.unique2[sel==1]
    #     batch.id.unique2 = batch.id.unique2[sel!=1]
    #     batch.tempj = batch.tempj[sel!=1]
    #   }
    #   batch.id.unique.new[n_batch] = batch.id.unique[!(batch.id.unique %in% batch.id.unique.new)]
    #   batch.obsj = batch.obs[,c('batch_id', g.names[j])] %>%
    #     arrange(batch_id = factor(batch_id, levels = batch.id.unique.new))
    #   batch.obs.sorted[,j] = batch.obsj[,2]
    # }
    
    ## sort batch by survival time
    if(beta_sort_train != 0) {
      temp = exp(beta_sort_train*log(t.train))
      id2 = 1:N.train
      id.new = rep(NA, N.train)
      
      for(i in 1:(N.train-1)){
        probs = temp/sum(temp)
        sel = rmultinom(1, 1, probs)
        id.new[i] = id2[sel==1]
        id2 = id2[sel!=1]
        temp = temp[sel!=1]
      }
      id.new[N.train] = id[!(id %in% id.new)]
      
      # Sort batch effect
      ### [Per Andy 7/21/2025] ----
      ### Should be “batch.train.log = batch.obs.sorted+ xtrue.train.log[id.new,]”,
      ### i.e. it is the xtrue.train.log that needs to be partially sorted by id.new, but batch.obs.sorted.
      ### This is because in line 618 you sort t.train by id.new. 
      batch.train.log = xtrue.train.log[id.new,] + batch.obs.sorted
      batch.train = exp(as.data.frame(lapply(batch.train.log, cap_and_add_noise, cap=rwd.max, sd=rwd.sd)))
    }
    
    # Final train data
    x.train.count = batch.train
    x.train = log(x.train.count+0.5)
  }

  ### HE Test -------------
  batch.id.test = batch.event.temp[test.ids,'batch.id']
  id.test = id.new.test = 1:N.test
  
  if(he_test==0) {
    x.test.count=xtrue.test
    x.test=log(x.test.count+0.5)
  
  } else {
    xtrue.test.log = log(xtrue.test+0.5)
    batch.test.log = log(batch.test+0.5)
    ## NOTE: on the natural log scale
    ### update 7/16: scale batch effects by 1/2
    batch.obs = (batch.test.log - xtrue.test.log) / 2 
    batch.test.log = batch.obs + xtrue.test.log
    batch.test = exp(as.data.frame(lapply(batch.test.log, cap_and_add_noise, cap=rwd.max, sd=rwd.sd)))
    
    ## per batch median
    batch.obs$batch_id = batch.id.test
    batch.order = cbind(
      Median = apply(batch.obs %>% group_by(batch_id) %>% summarise(across(everything(), median)) %>% select(-batch_id), 1, median),
      batch_id = seq(1,10)
    ) %>% as.data.frame() %>% 
      arrange(desc(Median)) %>% pull(batch_id)
    
    batch.obs$batch_id = factor(batch.obs$batch_id, levels = batch.order)
    batch.obs.sorted = batch.obs[order(batch.obs$batch_id), ] %>% select(-batch_id)
    
    # batch.obs.med = batch.obs %>%
    #   group_by(batch_id) %>%
    #   summarise(across(everything(), median)) %>%
    #   select(-batch_id) %>%
    #   as.matrix()
    # batch.temp = exp(batch_sort_test*batch.obs.med)
    # 
    # ## partially sort batch.obs for each gene
    # batch.obs.sorted = matrix(nrow=N.test, ncol=P)
    # colnames(batch.obs.sorted) = g.names
    # 
    # for(j in 1:P){
    #   batch.id.unique2=batch.id.unique
    #   batch.id.unique.new=rep(NA, n_batch)
    #   batch.tempj=batch.temp[,j]
    #   
    #   for(b in 1:(n_batch-1)){
    #     probs=batch.tempj/sum(batch.tempj)
    #     sel=rmultinom(1, 1, probs)
    #     batch.id.unique.new[b]=batch.id.unique2[sel==1]
    #     batch.id.unique2=batch.id.unique2[sel!=1]
    #     batch.tempj=batch.tempj[sel!=1]
    #   }
    #   batch.id.unique.new[n_batch]=batch.id.unique[!(batch.id.unique %in% batch.id.unique.new)]
    #   batch.obsj = batch.obs[,c('batch_id', g.names[j])] %>%
    #     arrange(batch_id=factor(batch_id, levels = batch.id.unique.new))
    #   batch.obs.sorted[,j]=batch.obsj[,2]
    # }

    ## Sort batch by survival time
    if(beta_sort_test!=0){
      temp=exp(beta_sort_test*log(t.test))
      id2=1:N.test
      id.new.test=rep(NA, N.test)
      
      for(i in 1:(N.test-1)) {
        probs = temp/sum(temp)
        sel = rmultinom(1, 1, probs)
        id.new.test[i] = id2[sel==1]
        id2 = id2[sel!=1]
        temp = temp[sel!=1]
      }
      id.new.test[N.test] = id.test[!(id.test %in% id.new.test)]
      
      # Sort batch effect
      batch.test.log = xtrue.test.log[id.new.test,] + batch.obs.sorted
      batch.test = exp(as.data.frame(lapply(batch.test.log, cap_and_add_noise, cap=rwd.max, sd=rwd.sd)))
    }
    
    # Final test data
    x.test.count = batch.test
    x.test = log(x.test.count+0.5)
  }

  
  ## 5. Normalization -------------------------------
  if(norm_type==1){ # total count
    print("Normalization type: TC")

    x.train_tc=norm.TC(raw=t(x.train.count))
    x.train.count=t(x.train_tc$dataNormalized)
    x.train=log(x.train.count+0.5)

    x.test_tc=norm.TC(raw=t(x.test.count))
    x.test.count=t(x.test_tc$dataNormalized)
    x.test=log(x.test.count+0.5)
  }

  if(norm_type==2) { # upper quantile
    print("Normalization type: Upper Quantile")

    x.train_uq=norm.UQ(raw=t(x.train.count))
    x.train.count=t(x.train_uq$dataNormalized)
    x.train=log(x.train.count+0.5)

    x.test_uq=norm.UQ(raw=t(x.test.count))
    x.test.count=t(x.test_uq$dataNormalized)
    x.test=log(x.test.count+0.5)
  }
  
  if(norm_type==3) { # tmm
    print("Normalization type: TMM")

    x.train_tmm=norm.TMM(raw=t(x.train.count))
    x.train.count=t(x.train_tmm$dataNormalized)
    x.train=log(x.train.count+0.5)

    x.test_tmm=norm.TMM(raw=t(x.test.count))
    x.test.count=t(x.test_tmm$dataNormalized)
    x.test=log(x.test.count+0.5)
  }

  if(norm_type==4) { # quantile normalization
    print("Normalization type: Quantile Normalization")

    x.train=t(normalize.quantiles(t(x.train)))
    x.train.count=exp(x.train)
    x_quantile=sort(x.train[1,])

    x.test=t(apply(x.test, 1, function(u){x_quantile[rank(u)]}))
    x.test.count=exp(x.test)
  }
  
  if(norm_type==5) { # DESeq2
    print("Normalization type: DESeq")

    x.train_deseq = norm.DESeq(round(t(x.train.count),0) + 1)
    x.train.count = t(x.train_deseq$dataNormalized)
    x.train = log(x.train.count+0.5)
    
    x.test_deseq = norm.DESeq(round(t(x.test.count),0) + 1)
    x.test.count = t(x.test_deseq$dataNormalized)
    x.test = log(x.test.count+0.5)
  }
  
  if(norm_type==6) { # median
    print("Normalization type: Median")
    
    x.train_med=norm.med(raw=t(x.train.count))
    x.train.count=t(x.train_med$dataNormalized)
    x.train=log(x.train.count+0.5)
    
    x.test_med=norm.med(raw=t(x.test.count))
    x.test.count=t(x.test_med$dataNormalized)
    x.test=log(x.test.count+0.5)
  }
  
  colnames(x.train)=colnames(x.train.count)=colnames(x.test)=colnames(x.test.count)=g.names
  
  
  ## 6. Save Survival Data --------------------------
  
  # Code scenario namee e.g., "BE11Asso00_normTC" ({batchType}_norm{normType})
  convert_num_to_indicator <- function(x) {
    case_when(
      x>0  ~ 1,
      x==0 ~ 0,
      x<0  ~ -1
    )
  }
  sort_train = convert_num_to_indicator(beta_sort_train)
  sort_test =  convert_num_to_indicator(beta_sort_test)
  batchNormType = glue::glue("BE{he_train}{he_test}Asso{sort_train}{sort_test}_norm{norm_map[norm_type+1]}")
  
  # Final train and test data
  sim.train = data.frame(batch.id=batch.id.train, time=t.train[id.new],     status=delta.train[id.new],     x.train)
  sim.test  = data.frame(batch.id=batch.id.test,  time=t.test[id.new.test], status=delta.test[id.new.test], x.test)
  out = list(sim.train, sim.test)
  
  out.dir = file.path('data', batchNormType, sim.dataType)
  dir.create(out.dir, showWarnings = F, recursive = T)
  
  if (save_gene_data) {
    # Save survival data
    write.csv(sim.train %>% select(-c(time, status)),
              file.path('data', batchNormType, glue::glue("simGeneExp_train_{N.train}_{date}.csv")),
              row.names = F)
    write.csv(sim.test %>% select(-c(time, status)),
              file.path('data', batchNormType, glue::glue("simGeneExp_test_{N.test}_{date}.csv")),
              row.names = F)
  }
  
  if (save_surv_data) {
    # Save gene expression w/ or w/out batch
    write.csv(sim.train %>% select(time, status),
              file.path('data', batchNormType, sim.dataType, 
                        glue::glue("simSurvival_{sim.dataType}_train_{N.train}_{date}.csv")),
              row.names = F)
    write.csv(sim.test %>% select(time, status),
              file.path('data', batchNormType, sim.dataType,
                        glue::glue("simSurvival_{sim.dataType}_test_{N.test}_{date}.csv")),
              row.names = F)
  }
  
  if (save_beta) {
    # Save true beta coefs
    write.csv(data.frame(beta0), 
              file.path(out.dir, glue::glue("beta0_{sim.dataType}_{date}.csv")) )
  }

  ## 7. K-M plot of simulated survival data -------
  
  if (plot_km) {
    sim.out=rbind(
      sim.train %>% dplyr::mutate(type='Train'),
      sim.test %>% dplyr::mutate(type='Test')
    )
    fit = survfit(Surv(time, status)~type, data=sim.out, conf.type="log-log")
    fig = ggsurvfit(fit) +
      labs(
        title = 'Survival',
        x = 'Time',
        y = 'Survival Probability'
      ) +
      coord_cartesian(xlim = c(0, max.censor.time), ylim = c(0,1)) +
      scale_x_continuous(breaks = seq(0, max.censor.time, by=48)) +
      scale_y_continuous(breaks = seq(0, 1, by=0.2)) +
      add_confidence_interval() +
      add_risktable() +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      theme_bw()
    
    out = list(sim.train, sim.test, fig)
  }
  
  
  ## 8. Analysis --------------------------------
  if (run_analysis) {
    
    cat("\nMoving on to CoxPH analysis data..\n\n")
    
    # Initialize metric vectors
    n_train=c_o=c_o_test=c_nl=c_nl_test=c_l=c_l_test=l_l=c()
    
    if(stratify == 1){
      c_os=c_os_test=c_nls=c_nls_test=c_ls=c_ls_test=l_ls=c()
    }
    
    # Convert train data to Python data frame for splitting
    sim_train_py = r_to_py(sim.train) 
    sim_train_py$strata = 
      sim_train_py$status$astype("str")$str$cat(
        sim_train_py[['batch.id']]$astype("str"), sep = "_"
      )
    
    # Convert test data to Python data frame
    sim_test_py = r_to_py(sim.test) 
    sim_test_py$strata = 
      sim_test_py$status$astype("str")$str$cat(
        sim_test_py[['batch.id']]$astype("str"), sep = "_"
      )
    
    for (i in 1:length(subset_sizes)) {
      # i = 2
      n = subset_sizes[i]
      n_run = runs_per_size[i]
      n_splits = splits_per_size[i]
      run_c_o=run_c_o_test=run_c_nl=run_c_nl_test=run_c_l=run_c_l_test=c()
      
      if (stratify == 1) {
        run_c_os=run_c_os_test=run_c_nls=run_c_nls_test=run_c_ls=run_c_ls_test=c()
      }
      
      cat(glue::glue("Running for N={n}...\n\n"))
      
      for (run in 1:n_run) {
      # run=1
        n_train[length(n_train)+1] = n
        
        # Sample training subset -------------------
        if (n < nrow(sim_train_py)) {
          train_sub_py = skl_ms$train_test_split(
            sim_train_py,
            train_size=as.integer(n),
            shuffle=TRUE,
            random_state=as.integer(run),
            stratify=sim_train_py$strata
          )[[1]]
          
        } else {
          train_sub_py = sim_train_py
        }
        # Sample testing subset -------------------------
        ## [update 6/3/2025] also vary test set across iterations
        test_sub_py = skl_ms$train_test_split(
          sim_test_py,
          train_size=as.integer(test_size),
          shuffle=TRUE,
          random_state=as.integer(run),
          stratify = sim_test_py$strata
        )[[1]]
        
        data_train = py_to_r(train_sub_py)
        data_test  = py_to_r(test_sub_py) #sim.test
        x_train = data_train[, !colnames(data_train) %in% c('batch.id','time','status')]
        x_test  = data_test[, !colnames(data_test) %in% c('batch.id','time','status')]
        x0_train = data_train[,selected.geneid]
        x0_test  = data_test[, selected.geneid]
        
        ### Non-stratified analysis ------------------------------------------
        
        #### Oracle (linear) analysis ------------------------------------
        # Train
        coxfit_o = tryCatch(
          coxph(Surv(data_train$time, data_train$status) ~ as.matrix(x0_train)),
          error=function(e) e,
          warning=function(w) w
        )
        if(is(coxfit_o, "warning") | is(coxfit_o, "error")) {
          coxfit_o = coxph(Surv(data_train$time, data_train$status) ~ ridge(as.matrix(x0_train), theta=1))
        }
        run_c_o[length(run_c_o)+1] = c_o[length(c_o)+1] =
          summary(coxfit_o)$concordance[1]

        # Test
        coxfit_o_test = coxph(
          Surv(data_test$time, data_test$status) ~ as.matrix(x0_test),
          init=summary(coxfit_o)$coefficient[,1],
          control=coxph.control(iter.max=0)
        )
        run_c_o_test[length(run_c_o_test)+1] = c_o_test[length(c_o_test)+1] =
          summary(coxfit_o_test)$concordance[1]
        
        #### Oracle (nonlinear) analysis --------------------------------
        # Transform train/test gene count data (rows: samples x columns: genes)
        if (!sim.dataType %in% c('linear-moderate', 'linear-weak')) {
          
          ##### 10 quadratic terms -------------------
          if (sim.dataType=='nl-quadratic') {
            
            x0_t_train = apply(x0_train, 2, scale)^2  ## scale across 2=columns (genes)
            x0_t_test = apply(x0_test, 2, scale)^2
            
          #### 10 quadratic with intercept ----------
          } else if (sim.dataType=='nl-shiftquad') {
            
            x0_train_scaled = apply(x0_train, 2, scale)
            x0_test_scaled = apply(x0_test, 2, scale)
            
            c_tr = apply(x0_train_scaled, 2, median)
            c_te = apply(x0_test_scaled, 2, median)
            
            x0_t_train = sweep(x0_train_scaled, 2, c_tr, "-")^2   # element-wise transformation
            x0_t_test = sweep(x0_test_scaled, 2, c_te, "-")^2
            
          ##### 10 interaction terms ----------------
          } else if (sim.dataType=='nl-interaction') {
         
            # cat(glue::glue("Interaction term: {interact.gene}\n"))
            x0_train_scaled = apply(x0_train, 2, scale)
            x0_test_scaled = apply(x0_test, 2, scale)
            interact.gene.tr = x0_train_scaled[,interact.gene]
            interact.gene.te = x0_test_scaled[, interact.gene]
            
            x0_t_train = sweep(x0_train_scaled, 1, interact.gene.tr, FUN="*")
            x0_t_test = sweep(x0_test_scaled, 1, interact.gene.te, FUN="*")
            
          ##### 10 sine terms ----------------------
          } else if (sim.dataType=='nl-sine') {
            
            x0_t_train = sin(x0_train) 
            x0_t_test = sin(x0_test)
          } 
          
          # Train
          coxfit_nl = tryCatch(
            coxph(Surv(data_train$time, data_train$status) ~ as.matrix(x0_t_train)), 
            error=function(e) e,
            warning=function(w) w
          )
          if(is(coxfit_nl, "warning") | is(coxfit_nl, "error")) {
            coxfit_nl = coxph(
              Surv(data_train$time, data_train$status) ~ ridge(as.matrix(x0_t_train), theta=1))
          }
          run_c_nl[length(run_c_nl)+1] = c_nl[length(c_nl)+1] = 
            summary(coxfit_nl)$concordance[1]      
          
          # Test
          coxfit_nl_test = coxph(
            Surv(data_test$time, data_test$status) ~ as.matrix(x0_t_test), 
            init=summary(coxfit_nl)$coefficient[,1], 
            control=coxph.control(iter.max=0)
          )
          run_c_nl_test[length(run_c_nl_test)+1] = c_nl_test[length(c_nl_test)+1] = 
            summary(coxfit_nl_test)$concordance[1]  
        
        } else {
          run_c_nl[length(run_c_nl)+1] = 
            c_nl[length(c_nl)+1] = 
            run_c_nl_test[length(run_c_nl_test)+1] =
            c_nl_test[length(c_nl_test)+1] = NA
        }   
        
        #### Penalized Lasso ---------------------------------------------
        # ## Univariate filtering
        # lkhd_u = rep(0, ncol(x_train))
        # for(j in 1:length(lkhd_u)){
        #   coxfit_ui=coxph(Surv(data_train$time, data_train$status) ~ as.matrix(x_train[,j]))
        #   if(!is.na(coef(coxfit_ui)) & abs(coef(coxfit_ui))<20){
        #     lkhd_u[j]=coxfit_ui$loglik[2]
        #   }else{
        #     lkhd_u[j]=-10000
        #   }
        # }
        # lkhd_p_sort=sort(lkhd_u, decreasing=TRUE)
        # lkhd_p_thres=lkhd_p_sort[round(min(sum(data_train$status)/10, ncol(x_train)/4))]
        
        cv_l = cv.glmnet(
          x_train,
          # as.matrix(x_train[,lkhd_u >= lkhd_p_thres]),
          Surv(data_train$time, data_train$status),
          nfolds=n_splits,
          family="cox",
          alpha=1, # Lasso
          standardize=F
          )
        lambda_min = cv_l$lambda.min
        l_l[length(l_l)+1] = lambda_min

        # Fit final lasso regression
        glmnet_l = glmnet(
          as.matrix(x_train),
          # as.matrix(x_train[,lkhd_u >= lkhd_p_thres]),
          Surv(data_train$time, data_train$status),
          family="cox", alpha=1, lambda=lambda_min,
          standardize=F
          )

        ## Train score
        lp_train = predict(glmnet_l, newx = as.matrix(x_train), #as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), 
                           s=lambda_min, type="link")
        c_l[length(c_l)+1] = run_c_l[length(run_c_l)+1] =
          survcomp::concordance.index(x=lp_train, surv.time=data_train$time, surv.event=data_train$status)$c.index

        ## Test score
        lp_test = predict(glmnet_l, newx = as.matrix(x_test), #as.matrix(x_test[,lkhd_u>=lkhd_p_thres]),
                          s=lambda_min, type="link")
        c_l_test[length(c_l_test)+1] = run_c_l_test[length(run_c_l_test)+1] =
          survcomp::concordance.index(x=lp_test, surv.time=data_test$time, surv.event=data_test$status)$c.index


        ### Stratified analysis ----------------------------------------------
        if(stratify == 1) {
          
          #### Oracle (linear) -------------------------
          # Train
          coxfit_os = tryCatch(
            coxph(Surv(data_train$time, data_train$status)~as.matrix(x0_train) + strata(data_train$batch.id)),
            error = function(e) e, 
            warning = function(w) w
          )
          if (is(coxfit_os, "warning") | is(coxfit_os, "error")){
            coxfit_os = 
              coxph(Surv(data_train$time, data_train$status)~ridge(as.matrix(x0_train), theta=1)+strata(data_train$batch.id))
          }
          run_c_os[length(run_c_os)+1] = c_os[length(c_os)+1] = summary(coxfit_os)$concordance[1]
          
          # Test
          coxfit_os_test = coxph(
            Surv(data_test$time, data_test$status)~as.matrix(x0_test) + strata(data_test$batch.id), 
            init = summary(coxfit_os)$coefficient[,1], 
            control = coxph.control(iter.max=0)
          )
          run_c_os_test[length(run_c_os_test)+1] = c_os_test[length(c_os_test)+1] =
            summary(coxfit_os_test)$concordance[1]

          #### Oracle (nonlinear) ------------------------
          if (!sim.dataType %in% c('linear-moderate', 'linear-weak')) {
            # Train
            coxfit_nls = tryCatch(
              coxph(Surv(data_train$time, data_train$status)~as.matrix(x0_t_train) + strata(data_train$batch.id)),
              error = function(e) e, 
              warning = function(w) w
            )
            if(is(coxfit_nls, "warning") | is(coxfit_nls, "error")){
              coxfit_os = 
                coxph(Surv(data_train$time, data_train$status)~ridge(as.matrix(x0_t_train), theta=1)+strata(data_train$batch.id))
            }
            run_c_nls[length(run_c_nls)+1] = c_nls[length(c_nls)+1] = summary(coxfit_nls)$concordance[1]      
            
            # Test
            coxfit_nls_test = coxph(
              Surv(data_test$time, data_test$status)~as.matrix(x0_t_test)+strata(data_test$batch.id), 
              init = summary(coxfit_nls)$coefficient[,1], 
              control = coxph.control(iter.max=0)
            )
            run_c_nls_test[length(run_c_nls_test)+1] = c_nls_test[length(c_nls_test)+1] = 
              summary(coxfit_nls_test)$concordance[1] 
          } else {
              run_c_nls[length(run_c_nls)+1] = 
                c_nls[length(c_nls)+1] = 
                run_c_nls_test[length(run_c_nls_test)+1] =
                c_nls_test[length(c_nls_test)+1] = NA
            }   
          }
          
          #### Lasso-penalized analysis with handling effect ------
          if (run > 5) {next}
        
          ## Univariate filtering
          lkhd_us = rep(0, ncol(x_train))
          for(j in 1:length(lkhd_us)){
            coxfit_usi=tryCatch(
              coxph(Surv(data_train$time, data_train$status) ~ as.matrix(x_train[,j])+strata(data_train$batch.id)), 
              error=function(e) e
              )
            if(!is(coxfit_usi, "error")) {
              if(!is.na(coef(coxfit_usi))) {
                if(abs(coef(coxfit_usi))<20) {
                  lkhd_us[j]=coxfit_usi$loglik[2]
                } else{
                  lkhd_us[j]=-10000
                }
              } else{
                lkhd_us[j]=-10000
              }
            } else{
              lkhd_us[j]=-10000
            }
          }
          lkhd_ps_sort=sort(lkhd_us, decreasing=TRUE)
          lkhd_ps_thres=lkhd_ps_sort[round(min(sum(data_train$status)/10, length(lkhd_us)/4))]
          
          ## Fit stratified Lasso CV regression 
          cv_ls = glmnet_strat_cv(
            data=cbind(
              data_train %>% tibble::rownames_to_column('id') %>% select(id,t=time,delta=status,batch.id),
              x_train[,lkhd_us>=lkhd_ps_thres]
              ), 
            lambda_max=0.3, 
            nlambda=2, 
            nfold=n_splits, 
            penalty_wt=rep(1, sum(lkhd_us>=lkhd_ps_thres)), 
            threshold=10^(-5)
            )
          b_ls0temp = cv_ls[[1]]
          l_ls[length(l_ls)+1] = cv_ls[[2]]
          
          # Train
          coxfit_ls = coxph(
            Surv(data_train$time, data_train$status)~as.matrix(x_train[,lkhd_us>=lkhd_ps_thres])+strata(data_train$batch.id), 
            init=b_ls0temp, 
            control=coxph.control(iter.max=0))
          
          run_c_ls[length(run_c_ls)+1]=c_ls[length(c_ls)+1]=
            summary(coxfit_ls)$concordance[1]
          
          # Test
          coxfit_ls_test = coxph(
            Surv(data_test$time, data_test$status) ~ as.matrix(x_test[,lkhd_us>=lkhd_ps_thres])+strata(data_test$batch.id), 
            init=b_ls0temp, 
            control=coxph.control(iter.max=0))
          
          run_c_ls_test[length(run_c_ls_test)+1]=c_ls_test[length(c_ls_test)+1]=
            summary(coxfit_ls_test)$concordance[1]
        } # end of multiple runs for loop
      
      cat(glue::glue("=============== Non-stratified ===============
    (Oracle linear)     Train: {round(mean(run_c_o,na.rm=T),3)} |  Test: {round(mean(run_c_o_test, na.rm=T),3)}\
    (Oracle nonlinaer)  Train: {round(mean(run_c_nl,na.rm=T),3)} |  Test: {round(mean(run_c_nl_test, na.rm=T),3)}\
    (Lasso)             Train: {round(mean(run_c_l, na.rm=T),3)} |  Test: {round(mean(run_c_l_test, na.rm=T),3)}\n\n"
      ))
      if (stratify == 1) {
        cat(glue::glue("================= Stratified ==================
      (Oracle linear)     Train: {round(mean(run_c_os,na.rm=T),3)} |  Test: {round(mean(run_c_os_test, na.rm=T),3)}\
      (Oracle nonlinaer)  Train: {round(mean(run_c_nls,na.rm=T),3)} |  Test: {round(mean(run_c_nls_test, na.rm=T),3)}\
      (Lasso)             Train: {round(mean(run_c_ls, na.rm=T),3)} |  Test: {round(mean(run_c_ls_test, na.rm=T),3)}\n\n"
        ))    
      }
    } # end of all subsets for loop
    
    ### Save model results ------------
    modelList = c('o','nl', 'l')
    if (stratify==1) { modelList = c(modelList, c('os','nls','ls')) }
    
    for(model in modelList) {
      
      if (model=='o') {
        modelName='oracle-linear'
      } else if (model=='nl') {
        modelName='oracle-nl'
      } else if (model=='l') {
        modelName='lasso'
      } else if (model=='os') {
        modelName='stratified-oracle-linear'
      } else if (model=='nl') {
        modelName='stratified-oracle-nl'
      } else {
        modelName='stratified-lasso'
      }
      
      results.dir = file.path("models", batchNormType, sim.dataType, modelName)
      dir.create(results.dir, showWarnings=F, recursive = T)
      
      cobj=eval(parse(text=paste0("c_",model)))
      cobj_test=eval(parse(text=paste0("c_",model,"_test")))
      model_results = data.frame('n train' = n_train, 
                                 'train C' = cobj, 
                                 'test C'  = cobj_test)
      write.csv(model_results, 
                file.path(results.dir, paste0("model_results_",Sys.Date(),".csv")), 
                row.names=F)
      
      if (model %in% c("l","ls")){
        lobj=eval(parse(text=paste0("l_", model)))
        write.csv(lobj,
                  file.path(results.dir, paste0("lambda_min_", Sys.Date(),".csv")),
                  row.names=F)      
      }
    }
  }
  
  # 9. Generate Slurm script for job submission -----------
  if (generate_subFile) {
    
    ## DL config file ----------
    
    config = list(
      batchNormType = batchNormType,
      dataName = sim.dataType,
      storage_url = "sqlite:///deepsurv-torch-hp-log-1.db",
      keywords = list(glue::glue("{date}")),
      time_col = "time",
      status_col = "status",
      batch_col = "batch.id",
      keep_batch = TRUE,
      hyperparameters = list(
        num_nodes = list(type = "categorical", choices = list(list(128), list(64), list(32), list(16),
                                                              c(128,128), c(64,64), c(32,32), c(16,16))),
        dropout = list(type = "float", low = 0.1, high = 0.5),
        learning_rate = list(type = "float", low = 1e-4, high = 5e-3, log = TRUE),
        weight_decay = list(type = "float", low = 1e-4, high = 5e-2, log = TRUE),
        batch_size = list(type = "categorical", choices = c(32, 64, 128))
      ),
      subset_sizes = c(100, 500, 1000, 2000, 5000, 10000),
      runs_per_size = c(20, 20, 20, 20, 20, 20),
      splits_per_size = c(3, 5, 5, 10, 10, 10),
      trials_per_size = c(25, 25, 25, 25, 25, 25),
      trial_threshold = 5,
      n_jobs = 1,
      is_tune = TRUE,
      is_save = TRUE,
      early_stop_per_size = list(
        `100` = list(patience = 20, min_delta = 5e-3),
        `500` = list(patience = 20, min_delta = 5e-3),
        `1000` = list(patience = 30, min_delta = 1e-3),
        `2000` = list(patience = 30, min_delta = 5e-4),
        `5000` = list(patience = 35, min_delta = 1e-4),
        `10000` = list(patience = 40, min_delta = 1e-4)
      )
    )
    # Save to json file
    file = here::here('configs', glue::glue("{batchNormType}-{sim.dataType}.json"))
    if (!file.exists(file)) {
      jsonlite::write_json(config, file, pretty=T, auto_unbox=T)
    }
               
    ## DL Slurm script --------------------
    
    bash.sub = glue::glue(
"#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}{GPUtype}
#SBATCH --error=jobs/deepsurv/logs/{batchNormType}-{sim.dataType}{GPUtype}.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition={tolower(GPUtype)}gpu
#SBATCH --qos={tolower(GPUtype)}gpu
#SBATCH --gres=gpu:{toupper(GPUtype)}:1
#SBATCH --time=02:00:00

export CUDA_HOME=/opt/cuda118
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$CUDA_HOMSE/bin:$PATH\n
source /home/nfs/dengy/dl-surv/bin/activate
python3 main.py --config configs/{batchNormType}-{sim.dataType}.json\n
deactivate"
    )
    file = here::here('jobs','deepsurv',glue::glue("{batchNormType}-{sim.dataType}.sh"))
    if (! file.exists(file)) {
      # dir.create(here::here('jobs','deepsurv'), showWarnings=F, recursive=T)
      write.table(bash.sub, file, row.names=F, col.names=F, quote=F)
    }

    ## ML bash script --------------------
    
    bash.sksurv = glue::glue(
"#!/bin/bash

LOGFILE='jobs/sksurv/logs/{batchNormType}-{sim.dataType}.log'
exec > >(tee -a $LOGFILE) 2>&1

source ~/dl-surv/bin/activate

COMMON_ARGS='--batchNormType {batchNormType} \\
--dataName {sim.dataType} \\
--subsets 100,500,1000,2000,5000,10000 \\
--runs 20,20,20,20,20,20 \\
--splits 3,5,5,10,10,10 \\
--trials 20,20,20,20,20,20 \\
--keywords {date} \\
--is_save'

python main-sksurv.py $COMMON_ARGS --modelType svm
python main-sksurv.py --batchNormType {batchNormType} \\
--dataName {sim.dataType} \\
--subsets 100,500,1000,2000,5000,10000 \\
--runs 20,20,20,20,5,5 \\
--keywords {date} \\
--is_save --modelType rsf
python main-sksurv.py $COMMON_ARGS --modelType gb

echo 'Congrats! All survival models completed.'"
    )
    file = here::here('jobs','sksurv',glue::glue("{batchNormType}-{sim.dataType}.sh"))
    
    if (!file.exists(file)) {
      # dir.create(here::here('jobs','sksurv'), showWarnings=F, recursive=T)
      write.table(bash.sksurv, file, row.names=F, col.names=F, quote=F)
    }
  }
  
  cat('Done.\n+++++++++++++++++++++++++++++++++++++++++++++++')
  
  return(out)
}
