library(ggsurvfit)
library(survival)
setwd("~/dl-survival-miRNA/")


# Load simulated miRNA -----
## Parameters ----
N = 5000
p = 10
date = format(Sys.Date(), "%m%d%y")
seed = 1234

########(subset N=5000 each) (log2+1 transformed)
mxf = read.csv(file.path('data',"MSKpair_MXF_5000.csv"))[1:N,1:1033]
pmfh = read.csv(file.path('data',"MSKpair_PMFH_5000.csv"))[1:N,1:1033]
g.names = read.csv("data/MSKpair_miRNA_names.csv")$mirna.names[1:1033] ## gene names
working=rbind(mxf, pmfh)
colnames(working)=g.names
# write.csv(working, "data/simulate_survival_10000.csv", row.names=F)


# Function for quick simulation -----------
simulate_T <- function(h0, h.log, n_train, max.censor.time) {
  # simulate uncensored survival time
  t0=-log(runif(n_train))*h0 / exp(h.log)
  
  # create censoring time
  ct=runif(n_train, 0, max.censor.time) ## generate censor time (capped by a max time)
  t=ifelse(t0<=ct, t0, ct)
  delta=1*(t0<=ct)
  ncase.train=sum(delta)
  
  hist(t) ## visualize (censored) survival time
  print(paste0("Event rate: ", ncase.train / n_train)) ## about 75% event rate
  summary(t) ## summarize simulated survival time
  return(list(t, delta))
}



sim.survdata <- function(surv.data, sim.dataType, sim.method, is_log_scale=T, 
                         cap.gene.value=3350000, h0=50, max.censor.time=200, 
                         beta0=NULL, simulate_beta=T, p=10, p.weak=30,
                         plot_km=T, save_data=T, file.name=NULL, save_beta=T,
                         seed=1234) {
  
  ## Prepare gene expression data -------------------------
  if (is_log_scale){
    workingt = t( exp(surv.data) ) # transform back to normal scale
  } else {
    workingt = t(surv.data)
  }
  workingt[workingt>cap.gene.value] = cap.gene.value + 
    rnorm(sum(workingt>cap.gene.value), 0, 100) # cap values greater than max in real sarcoma data
  
  
  ## Simulate Beta0's from augmented gene data ------------
  if (simulate_beta) {
    ## Filter genes by count per million
    cpm=cpm(workingt)
    keep=(rowSums(cpm>2)>=5)
    workingt=workingt[keep,]
    
    ## Select p true genes associated with survival time
    tc=rowSums(workingt)
    mid.tc=(tc>median(tc) & tc<quantile(tc, 0.75)) # select genes with total count in the third quartile
    cv=apply(workingt, 1, sd)/rowMeans(workingt)
    mid.cv=(cv > quantile(cv, 0.75)) # select genes with coefficient of variation (SD/mean) in the upper half
    geneid=row.names(workingt)
    candidate.geneid=geneid[mid.tc & mid.cv]
    
    # Select p (defualt=10) true genes
    set.seed(seed) # for reproducible results
    if (p<=length(candidate.geneid)) 
    {
      selected.geneid=sample(candidate.geneid, p, replace=F)
    } else 
    {
      selected.geneid = candidate.geneid
    }
  }
  
  ### Alternatively import true gene and beta0 coefficients ----
  else {
    selected.geneid = rownames(beta0)
  }
  
  
  ## Simulate survival outcome ---------------------------
  
  # log-transform selected true genes
  xtrue = workingt[selected.geneid,]
  xtrue.log = log2(xtrue+0.5)
  n_train = ncol(xtrue)
  
  
  ### Linear Simulation ------------------
  if (sim.dataType == 'linear') 
  {
    #### Moderate signals -------
    if (sim.method=='moderate') {
      
      print("Simulating from linear risk function with moderate signals")
      
      ## calculate beta's as inverse std of gene expression
      sd.genes = apply(xtrue.log, 1, sd)
      beta0 = rep(c(1,-1),p/2) * round(1/sd.genes, 1)
      
      ## calculate log risk
      h.log= t(xtrue.log) %*% beta0 
    }
    
    ### Weak Signals ------------
    else if (sim.method=='weak') {
      
      print("Simulating from linear risk function with weak signals")
      
      ## re-select (default=30) true genes associated with survival time
      set.seed(seed)
      if (p.weak <= length(candidate.geneid)) {
        selected.weak.geneid = sample(candidate.geneid, p.weak, replace=F)
      } else {
        selected.weak.geneid = candidate.geneid
      }
      
      # Log-transform gene expression
      xtrue.weak = workingt[selected.weak.geneid,]
      xtrue.weak.log = log2(xtrue.weak+0.5)
      
      # Simulate beta values
      sd.genes.weak = apply(xtrue.weak.log, 1, sd)
      #### Note: shrink beta coefficients to conrol for total signal
      beta.scale = p/p.weak
      beta0 = rep(c(1,-1),p/2) * round(1/sd.genes.weak*beta.scale, 1)
      
      ## calculate log risk
      h.log= t(xtrue.weak.log) %*% beta0
    }
    
    surv.out=simulate_T(h0, h.log, n_train, max.censor.time)
  }
  
  
  ### Nonlinear Simulation ----------------
  else if (sim.dataType=='nonlinear') {
    
    #### 10 quadratic terms --------
    if (sim.method=='nl-quadratic') {
      
      print("Simulating with quadratic transformation")
      
      ## Scale gene expr data
      xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
      xtrue.nl.log = as.matrix(xtrue.log.scaled)^2
      
      sd.genes = apply(xtrue.nl.log, 1, sd)
      beta0 = rep(c(1,-1),5)*round(5/sd.genes, 1)
      
      h.log = t(xtrue.nl.log) %*% beta0
    }
    
    
    ### 10 quadratic interaction terms ------
    else if (sim.method=='nl-interaction') {
      
      print("Simulating with non-linear interaction terms")
      
      set.seed(seed) # randomly select a interaction gene
      interact.gene = sample(selected.geneid, 1)
      print(paste0("Interaction term selected: ", interact.gene) )
      
      xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
      interact.gene.val = xtrue.log.scaled[interact.gene,]
      
      xtrue.nl.log = as.matrix(xtrue.log.scaled) %*% diag(interact.gene.val)
      sd.genes = apply(xtrue.nl.log, 1, sd)
      beta0 = rep(c(1,-1),5) * round(3/sd.genes, 1) 
      
      h.log = t(xtrue.nl.log) %*% beta0
    }
    
    
    ### 10 sine terms -----------
    else if (sim.method=='nl-sine') {
      
      print("Simulating with non-linear Sine interaction terms")
      
      xtrue.nl.log = t(sin(t(xtrue.log))) # sine transformation
      
      # calculate beta coefficients
      sd.genes = apply(xtrue.nl.log, 1, sd)
      beta0 = rep(c(1,-1),5)*round(1/sd.genes/2, 2)
      
      # calculate log risk
      h.log = t(xtrue.nl.log) %*% beta0
    }
    
    ### Gaussian log risk function -----
    # else if (sim.method=='nl-gauss') 
    # {
    #   print("Simulating with non-linear Gaussian log-risk function")
    #   xtrue.log.scaled = t(scale(t(xtrue.log)))
    #   x.mean = apply(xtrue.log.scaled, 1, mean)
    #   sigma  = 1.0
    #   h.log  = apply(xtrue.log.scaled, 2,
    #                  function(x) exp( -sum(abs(x-x.mean)^2) / (2*sigma^2) )
    #   )
    # }
    
    ## Save simulated survival time and censoring outcome
    surv.out = simulate_T(h0, h.log, n_train, max.censor.time)
    
  }
  
  
  ## Save survival data ---------------------------
  t = surv.out[[1]]
  delta = surv.out[[2]]
  
  sim.survival.out = data.frame(time=t, status=delta)
  colnames(sim.survival.out) <- c('time','status')
  sim.out = cbind(surv.data, sim.survival.out)
  out=sim.survival.out
  
  if (save_data) {
    out.dir  = file.path('data', sim.dataType)
    dir.create(out.dir, showWarnings = F)
    if (is.null(file.name)) {
      file.name = file.path(paste0("simSurvival_",sim.method,"_",N*2,"_",date,".csv"))
    }
    ## save simulated survival outcomes
    write.csv(out, file=file.path(out.dir, file.name), row.names=F)
    
    ## save beta0's
    if (save_beta) {
      write.csv(data.frame(beta0),
                file.path(out.dir, paste0("beta0_",sim.method,"_", date,".csv")) )
    }
  }
  
  
  
  ## K-M plot of simulated survival data -------
  if (plot_km) {
    
    fit = survfit(Surv(time, status)~1, data=sim.out, conf.type="log-log")
    fig <- ggsurvfit(fit) +
      labs(
        title = 'Survival',
        x = 'Time (Days)',
        y = 'Survival Probability'
      ) +
      coord_cartesian(xlim = c(0, max.censor.time), ylim = c(0,1)) +
      scale_x_continuous(breaks = seq(0, max.censor.time, by=20)) +
      scale_y_continuous(breaks = seq(0, 1, by=0.2)) +
      add_confidence_interval() +
      add_risktable() +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      theme_bw()
    fig
    out = list(sim.survival.out, fig)
  }
  
  return(out)
}



# Run Simulation ---------------------------------

# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data=working, 
                          sim.dataType='linear', sim.method='moderate',
                          h0=39, save_data = T, save_beta = T, seed=1234,
                          file.name=paste0("simSurvival_linear-moderate_10000_latest_RW.scale.csv"))
linear.moderate.surv = linear.moderate.out[[1]]
## Look at the distribution of simulated data
plot(survfit(Surv(linear.moderate.surv$time, linear.moderate.surv$status)~1))


# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data=working, 
                       sim.dataType='linear', sim.method='weak',
                       h0=40, save_data = T, save_beta = T, seed=1234,
                       file.name=paste0("simSurvival_linear-weak_10000_latest_RW.scale.csv"))
linear.weak.surv = linear.weak.out[[1]]


# Squared terms
nl.qua.out = sim.survdata(surv.data = working,
                        sim.dataType='nonlinear', sim.method='nl-quadratic',
                        h0=1.75, save_data = T, save_beta = T, seed = 1234,
                        file.name=paste0("simSurvival_nl-quadratic_10000_latest_RW.scale.csv"))
nl.qua.surv=nl.qua.out[[1]]


# Quadratic interaction
nl.interact.out = sim.survdata(surv.data=working, 
                        sim.dataType='nonlinear', sim.method='nl-interaction',
                        h0=49, save_data = T, save_beta = T, seed = 1234,
                        file.name=paste0("simSurvival_nl-interaction_10000_latest_RW.scale.csv"))
nl.interact.surv = nl.interact.out[[1]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data=working,
                          sim.dataType='nonlinear', sim.method='nl-sine',
                          h0=8.75, save_data = T, save_beta = T, seed = 1234,
                          file.name=paste0("simSurvival_nl-sine_10000_latest_RW.scale.csv"))
nl.sine.surv = nl.sine.out[[1]]


# # Gaussian kernel
# nl.gauss.rw.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-gauss',
#                              h0=140, simulate_beta = F, save_data = T,
#                              file.name=paste0("simSurvival_nl-gauss_10000_", date, "_RW.scale.csv"))
# nl.gauss.rw.surv = nl.gauss.rw.out[[1]]




# # [SKIP] Run Simulation with unscaled betas -------
# linear.exp.out = sim.survdata(surv.data=working, sim.dataType = 'linear', sim.method = 'exponential',
#                                  h0=0.4e-3, beta0=beta0, save_data = F, seed=1234,
#                                  file.name="simSurvival_exponential_10000_091824_orig.scale.csv")
# linear.exp.surv= linear.exp.out[[1]]
# 
# 
# nl.1.out = sim.survdata(surv.data = working, sim.dataType = 'nonlinear', sim.method = 'nl-1',
#                            h0=1.7e-2, beta0=beta0, save_data = F, seed=1234,
#                            file.name="simSurvival_nl-1_10000_091824_orig.scale.csv")
# nl.1.surv=nl.1.out[[1]]
# 
# 
# nl.qua.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method = 'nl-quadratic',
#                              h0=4e-5, beta0=beta0, save_data = F,
#                              file.name="simSurvival_nl-quadratic_10000_091824_orig.scale.csv")
# nl.qua.surv = nl.qua.out[[1]]
# 
# 
# nl.sine.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-sine',
#                               h0=1.5e-4, beta0=beta0, save_data = F,
#                               file.name="simSurvival_nl-sine_10000_091824_orig.scale.csv")
# nl.sine.surv = nl.sine.out[[1]]
