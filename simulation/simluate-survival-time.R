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
                         cap.gene.value=3350000, h0=150, max.censor.time=200, 
                         beta0=NULL, simulate_beta=F, p=NULL, save_beta_dir=NULL,
                         plot_km=T, save_data=T, file.name=NULL, seed=1234) {
  
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
    mid.cv=(cv>median(cv)) # select genes with coefficient of variation (SD/mean) in the upper half
    geneid=row.names(workingt)
    candidate.geneid=geneid[mid.tc & mid.cv]
    
    set.seed(seed) # for reproducible results
    if (p<=length(candidate.geneid)) 
    {
      selected.geneid=sample(candidate.geneid, p, replace=F)
    } else 
    {
      selected.geneid = candidate.geneid
    }
    
    ## Calculate beta's 
    selected.gene.cpm=cpm[selected.geneid,]
    sd.genes.cpm = apply(log2(selected.gene.cpm+0.5), 1, sd)
    
    ## Scale stdev of gene expression by constant c to obtain beta 0's
    #####Note: Maybe first log transform stdev's to shrink range (per Andriy)
    beta0=rep(c(1,-1),p/2)*round(1/sd.genes.cpm, 1)
    
    if (!is.null(save_beta_dir))
    {
      write.csv(data.frame(beta0),
        file.path('data', sim.dataType, paste0("beta0_",sim.method,"_", date,".csv")) )
    }
  }
  
  
  ## Simulate survival outcome -----------------------
  
  # Log-transform gene expression
  selected.geneid = rownames(beta0)
  xtrue = workingt[selected.geneid,]
  xtrue.log = log2(xtrue+0.5)
  n_train = ncol(xtrue)
  
  ### Linear Simulation ----------------
  if (sim.dataType == 'linear') {

    if (sim.method=='exponential') 
    {
      print("Simulating from exponential distribution")
      h.log= t(xtrue.log) %*% beta0 
      surv.out=simulate_T(h0, h.log, n_train, max.censor.time)
    }
  }
  ### Nonlinear Simulation -------------
  else if (sim.dataType=='nonlinear') {
  
    #### 1 quadratic term -------
    if (sim.method=='nl-1') 
    {
      print("Simulating with 1 quadratic interaction terms")

      set.seed(seed) # randomly select an interaction gene
      interact.gene = sample(selected.geneid, 1)
      print(paste0("Interaction term selected: ", interact.gene))# "hsa.miR.887.1"
      
      # [10/04] standardize gene expr before taking quadratic
      interact.col = (scale(xtrue.log[interact.gene,]))^2
      interact.beta=round(1/sd(interact.col), 1)
      
      xtrue.nl.log = rbind(xtrue.log, t(interact.col))
      beta0.nl = c(beta0.scaled, interact.beta)

      h.log = t(xtrue.nl.log) %*% beta0.nl
    }
    
    ### 10 quadratic terms --------
    else if (sim.method=='nl-quadratic')
    {
      print("Simulating with 10 quadratic interaction terms")
      
      set.seed(seed) # randomly select a interaction gene
      interact.gene = sample(selected.geneid, 1)
      print(paste0("Interaction term selected: ", interact.gene) )
      
      ######### [10/04] rescale X gene expression (std=1)
      gene.val.scaled = apply(xtrue.log[selected.geneid,], 1, scale)
      
      interact.col = as.matrix(gene.val.scaled)*as.vector(
                                  scale(xtrue.log[interact.gene,], center=F))
      interact.col = t(interact.col)
      
      sd.interact  = apply(interact.col, 1, sd)
      interact.beta= rep(c(1,-1),length(selected.geneid)/2)*round(1/sd.interact, 1)
      xtrue.nl.log = rbind(xtrue.log, interact.col)
      beta0.nl     = c(beta0.scaled, interact.beta)
      
      h.log = t(xtrue.nl.log) %*% beta0.nl
    }
    
    ### 10 sine terms -----------
    else if (sim.method=='nl-sine')
    {
      print("Simulating with non-linear Sine interaction terms")
      ####[10/04] drop randomness in sine simulation
      # err = matrix(rnorm(nrow(xtrue)*ncol(xtrue), 0, 0.25), ncol=nrow(xtrue)) 
      interact.col = sin(t(xtrue.log)) #+ err 

      sd.interact = apply(interact.col ,2, sd)
      interact.beta= rep(c(1,-1), length(selected.geneid)/2)*round(1/sd.interact, 1)
    
      xtrue.nl.log = rbind(xtrue.log, t(interact.col))
      beta0.nl = c(beta0, interact.beta)
      
      h.log = t(xtrue.nl.log) %*% beta0.nl
    }
    
    ### Gaussian log risk function -----
    else if (sim.method=='nl-gauss') 
    {
      print("Simulating with non-linear Gaussian log-risk function")
      xtrue.log.scaled = t(scale(t(xtrue.log)))
      x.mean = apply(xtrue.log.scaled, 1, mean)
      gamma  = 1.0
      h.log  = apply(xtrue.log.scaled, 2,
                    function(x) exp( -sum((x-x.mean)^2) * gamma )
      )
    }
    
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
    if (is.null(file.name)) 
    {
      file.name = file.path(paste0("simSurvival_",sim.method,"_",N*2,"_",date,".csv"))
    }
    write.csv(out, file=file.path(out.dir, file.name), row.names=F) ## save simulated survival outcomes
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

date
# Run Simulation -------------
# Linear exponential risk
linear.exp.rw.out = sim.survdata(surv.data=working, sim.dataType = 'linear', sim.method = 'exponential',
                        h0=23, beta0=beta0.scaled , save_data = T, seed=1234,
                        file.name=paste0("simSurvival_exponential_10000_", date, "_RW.scale.csv"))
linear.exp.rw.surv= linear.exp.rw.out[[1]]
# ## Look at the distribution of simulated data
# plot(survfit(Surv(linear.exp.rw.surv$time,linear.exp.rw.surv$status)~1))


# 1 quadratic interaction
nl.1.rw.out = sim.survdata(surv.data = working, sim.dataType = 'nonlinear', sim.method = 'nl-1',
                        h0=35, beta0=beta0.scaled , save_data = T,
                        file.name=paste0("simSurvival_nl-1_10000_", date, "_RW.scale.csv"))
nl.1.rw.surv=nl.1.rw.out[[1]]

# Quadratic interaction
nl.qua.rw.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method = 'nl-quadratic',
                           h0=8, beta0=beta0.scaled, save_data = T,
                           file.name=paste0("simSurvival_nl-quadratic_10000_", date, "_RW.scale.csv"))
nl.qua.rw.surv = nl.qua.rw.out[[1]]

# Sine interaction
nl.sine.rw.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-sine',
                          h0=3.8, beta0=beta0.scaled, save_data = T,
                          file.name=paste0("simSurvival_nl-sine_10000_", date, "_RW.scale.csv"))
nl.sine.rw.surv = nl.sine.rw.out[[1]]

# Gaussian kernel
nl.gauss.rw.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-gauss',
                             h0=140, simulate_beta = F, save_data = T,
                             file.name=paste0("simSurvival_nl-gauss_10000_", date, "_RW.scale.csv"))
nl.gauss.rw.surv = nl.gauss.rw.out[[1]]




# Run Simulation with unscaled betas -----
linear.exp.out = sim.survdata(surv.data=working, sim.dataType = 'linear', sim.method = 'exponential',
                                 h0=0.4e-3, beta0=beta0, save_data = F, seed=1234,
                                 file.name="simSurvival_exponential_10000_091824_orig.scale.csv")
linear.exp.surv= linear.exp.out[[1]]


nl.1.out = sim.survdata(surv.data = working, sim.dataType = 'nonlinear', sim.method = 'nl-1',
                           h0=1.7e-2, beta0=beta0, save_data = F, seed=1234,
                           file.name="simSurvival_nl-1_10000_091824_orig.scale.csv")
nl.1.surv=nl.1.out[[1]]


nl.qua.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method = 'nl-quadratic',
                             h0=4e-5, beta0=beta0, save_data = F,
                             file.name="simSurvival_nl-quadratic_10000_091824_orig.scale.csv")
nl.qua.surv = nl.qua.out[[1]]


nl.sine.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-sine',
                              h0=1.5e-4, beta0=beta0, save_data = F,
                              file.name="simSurvival_nl-sine_10000_091824_orig.scale.csv")
nl.sine.surv = nl.sine.out[[1]]
