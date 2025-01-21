library(edgeR)
library(ggsurvfit)
setwd("~/dl-survival-miRNA/")


# Load simulated miRNA -----
## Parameters
N = 5000
date = format(Sys.Date(), "%m%d%y")
seed = 1234

######## (subset N=5000 each) (log2+1 transformed)
# mxf = read.csv(file.path('data',"MSKpair_MXF_5000.csv"))[1:N,1:1033]
# pmfh = read.csv(file.path('data',"MSKpair_PMFH_5000.csv"))[1:N,1:1033]
mxf = read.csv(file.path('data',"MSKpair_MXF_AEhead_maf_60000_Draw1.csv"))[1:N,1:1033]
pmfh = read.csv(file.path('data',"MSKpair_PMFH_AEhead_maf_60000_Draw1.csv"))[1:N,1:1033]
g.names = read.csv("data/MSKpair_miRNA_names.csv")$mirna.names[1:1033] ## gene names
working=rbind(mxf, pmfh)
colnames(working)=g.names
# write.csv(working, "data/simulate_survival_10000.csv",row.names=F)


# ----- TO COLLAPSE ----------
workingt=t( exp(working) ) # transform back to normal scale
workingt[workingt>3350000] = 3350000 + 
  rnorm(sum(workingt>3350000),0,100) # cap values greater than max in real sarcoma data

### Simulate Beta0's from augmented gene data -----
## Filter genes by count per million
cpm=cpm(workingt)
keep=(rowSums(cpm>2)>=10) # results in 645 genes
# hist(rowSums(cpm>2), breaks=20)
workingt=workingt[keep,]


## Select p true genes associated with survival time
tc=rowSums(workingt)
mid.tc=(tc>median(tc) & tc<quantile(tc, 0.75)) # select genes with total count in the third quartile
cv=apply(workingt, 1, sd)/rowMeans(workingt)
# mid.cv=(cv>median(cv)) # select genes with coefficient of variation (SD/mean) in the upper half
mid.cv = (cv>quantile(cv, 0.75))
geneid=row.names(workingt)
candidate.geneid=geneid[mid.tc & mid.cv] # results in 73 genes


### Directly use Beta0's from RWD -----
beta0.scaled = beta0 = as.matrix(read.csv("data/beta0.csv", row.names = 'X'))
# "hsa.miR.1277.3p.1",0.2
# "hsa.miR.1277.5p.1",-0.3
# "hsa.miR.133a.2..1",0.2
# "hsa.miR.144..1",-0.2
# "hsa.miR.148b..1",0.5
# "hsa.miR.181c..1",-0.2
# "hsa.miR.204.1",0.2
# "hsa.miR.33a.1",-0.4
# "hsa.miR.598.1",0.3
# "hsa.miR.887.1",-0.2
# beta0[,1] = c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)



# Simulate survival outcome -----

set.seed(seed) # for reproducible results
p = 10
if (p<=length(candidate.geneid)) {
  selected.geneid=sample(candidate.geneid, p, replace=F)
} else {
  selected.geneid = candidate.geneid
}
selected.geneid
# "hsa.miR.301a..1"     "hsa.miR.525.5p.1"    "hsa.miR.190b..1"     "hsa.miR.892a.1"     
# "hsa.miR.129.5p.2..1" "hsa.miR.124.3..1"    "hsa.miR.329..2..1"   "hsa.miR.137.1"      
# "hsa.miR.1185.2.3p.1" "hsa.miR.599.1"  
# # New criteria:
# [1] "hsa.miR.3140.3p.1" "hsa.miR.934.1"     "hsa.miR.200b..1"   "hsa.miR.1293.3p.1"
# [5] "hsa.miR.1258.3p.1" "hsa.miR.370.1"     "hsa.miR.138.2..3"  "hsa.miR.124.3..1" 
# [9] "hsa.miR.659..1"    "hsa.miR.135b..1" 

## Log-transform gene expression
xtrue = workingt[selected.geneid,]
xtrue.log = log2(xtrue+0.5)
n_train = ncol(xtrue)


### TO COLLAPSE: Simsurv pacakge ----
# install.packages('simsurv');library(simsurv)
# # Example
# x=data.frame(t(xtrue.log))
# betas = data.frame(t(beta0.scaled))
# set.seed(seed)
# t <- simsurv(dist='exponential', lambdas=1/h0, x=x, #data.frame(hsa.miR.1277.3p.1=x$hsa.miR.1277.3p.1),
#              betas = c(
#                `hsa.miR.1277.3p.1`=0.2,
#                `hsa.miR.1277.5p.1`=-0.3,
#                `hsa.miR.133a.2..1`=0.2,
#                `hsa.miR.144..1`=-0.2,
#                `hsa.miR.148b..1`=0.5,
#                `hsa.miR.181c..1`=-0.2,
#                `hsa.miR.204.1`=0.2,
#                `hsa.miR.33a.1`=-0.4,
#                `hsa.miR.598.1`=0.3,
#                `hsa.miR.887.1`=-0.2)
# )
# t0=t$eventtime
# hist(t0)
# t$eventtime=ifelse(t0<=ct, t0, ct)
# t$status=1*(t0<=ct)  ## about 75% event rate
# t |> count(status)
# hist(t$eventtime)
# summary(t$eventtime)  
# 
# t.exp.ranked = t |> arrange(eventtime) |> tibble::rownames_to_column(var='rank') |> 
#   select(id, rank_simsurv=rank, eventtime_simsurv=eventtime, status_simsurv=status)
# t.exp.ranked.combined = left_join(t.exp.ranked, t.exp.rw.surv, by='id')

#### Cross reference 
# t.exp.surv = read.csv("data/linear/simulate_survival_10000_082924.csv")[,c('time','status')]
# t.exp.surv$id = seq(1, n_train)
# t.exp.surv =  t.exp.surv|>
#   arrange(time) |>tibble::rownames_to_column(var='rank')# |> tibble::column_to_rownames(var='idx')
# View(linear.exp.surv.ranked)
# 
# t.exp.ranked = left_join(t.exp.ranked.combined, t.exp.surv, by='id')
# # select(idx, time, status, rank, time_rw, status_rw, rank_rw)
# View(linear.exp.combined)




## Linear simulation -----
dataType = 'linear'

### Exponential: moderate =====
## Calculate beta's 
sd.genes = apply(xtrue.log, 1, sd)
sd.genes

beta0.moderate=rep(c(1,-1),p/2)*round(1/sd.genes, 1)
beta0.moderate

###### Note: modified to skip drawing from exponential distributions and directly
###### use the calculated hazard (h) as survival time (non-censored)
set.seed(seed)
h0=23
t0=-log(runif(n_train)) *h0 / exp(t(xtrue.mod.log)%*%beta0.moderate)
hist(t0[t0<200])

ct=runif(n_train, 0, 200) ## generate censoring time from unif(0,200)
t=ifelse(t0<=ct, t0, ct)
delta=1*(t0<=ct)  ## about 75% event rate
ncase.train=sum(delta)
ncase.train / n_train ## event rate
summary(t) ## simulated survival time
file.name = "simSurvival_exponential_moderate_10000_latest_RW.scale.csv"


### Exponential: Weak signal =====
method='weak'
set.seed(seed) # for reproducible results
p = 30
if (p<=length(candidate.geneid)) {
  selected.geneid=sample(candidate.geneid, p, replace=F)
} else {
  selected.geneid = candidate.geneid
}
# [1] "hsa.miR.301a..1"     "hsa.miR.525.5p.1"    "hsa.miR.190b..1"     "hsa.miR.892a.1"      "hsa.miR.129.5p.2..1" "hsa.miR.124.3..1"   
# [7] "hsa.miR.329..2..1"   "hsa.miR.137.1"       "hsa.miR.1185.2.3p.1" "hsa.miR.599.1"       "hsa.miR.659..1"      "hsa.miR.514a..3..1" 
# [13] "hsa.miR.520h.1"      "hsa.miR.520d.5p.1"   "hsa.miR.133b.1"      "hsa.miR.449b.1"      "hsa.miR.494.1"       "hsa.miR.890.3p.1"   
# [19] "hsa.miR.655..1"      "hsa.miR.182..1"      "hsa.miR.337.1"       "hsa.miR.551b..1"     "hsa.miR.665.1"       "hsa.miR.509.5p.2..1"
# [25] "hsa.miR.934.1"       "hsa.miR.506.1"       "hsa.miR.376a.2.5p.1" "hsa.miR.889.1"       "hsa.miR.376b..1"     "hsa.miR.1185.1.3p.1"

# Log-transform gene expression
xtrue.weak = workingt[selected.geneid,]
xtrue.weak.log = log2(xtrue.weak+0.5)

# Simulate beta values
sd.genes.weak = apply(xtrue.weak.log, 1, sd)
## Scale stdev of gene expression by constant c to obtain beta 0's
#####Note: first log transform stdev's to shrink range (per Andriy)
beta0.weak=rep(c(1,-1),p/2)*round(1/sd.genes.weak/3, 1)
beta0.weak
# write.csv(beta0.weak, file.path("data","beta0.weak.csv"), quote = F)
sum(abs(beta0.weak))
sum(abs(beta0.moderate))

h0=66
t0=-log(runif(n_train)) *h0 / exp(t(xtrue.weak.log)%*%beta0.weak)
hist(t0[t0<200])

ct=runif(n_train, 0, 200) ## generate censoring time from unif(0,200)
t=ifelse(t0<=ct, t0, ct)
delta=1*(t0<=ct)  ## about 75% event rate
ncase.train=sum(delta)
ncase.train / n_train ## event rate
summary(t) ## simulated survival time
file.name = "simSurvival_exponential_weak_10000_latest_RW.scale.csv"



## Nonlinear simulation -----
dataType='nonlinear'

#### non-linear 10 terms -----
method = "nl-quadratic"

## Scale gene expr data
xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
xtrue.nl2.log = as.matrix(xtrue.log.scaled)^2
View(xtrue.nl2.log)

sd.genes = apply(xtrue.nl2.log, 1, sd)
beta0.nl2= rep(c(1,-1),5)*round(5/sd.genes, 1)
beta0.nl2

h0=40
t0=-log(runif(n_train)) * h0*exp(t(interact.col)%*%interact.beta)
hist(t0[t0<200])
ct=runif(n_train, 0, 200) ## generate censoring time from unif(0,200)
t=ifelse(t0<=ct, t0, ct)
delta=1*(t0<=ct)  ## about 75% event rate
ncase.train=sum(delta)
print(paste0("Event rate: ", ncase.train / n_train))
summary(t)## simulated survival time 
file.name = "simSurvival_nl-quadratic_10000_latest_RW.scale.csv"


#### 10 non-linear terms: -----
method="nl-interaction"

# scale gene expression data
xtrue.log.scaled = t(apply(xtrue.log, 1, scale))

n.interact = 1
set.seed(seed) # randomly select a interaction gene
interact.gene = sample(selected.geneid, n.interact) #"hsa.miR.135b..1"
interact.gene.val = xtrue.log.scaled[interact.gene,]

xtrue.nl1.log = as.matrix(xtrue.log.scaled) %*% diag(interact.gene.val)
View(xtrue.nl1.log)
sd.interact = apply(xtrue.nl1.log, 1, sd)
# sd.interact
beta0.nl1 = rep(c(1,-1),5)*round(3/sd.interact, 1) 
beta0.nl1


h0=.7e-2
t0=-log(runif(n_train))*h0*exp(t(xtrue.nl1.log)%*%beta0.nl1)
hist(t0[t0<200])
ct=runif(n_train, 0, 200) ## generate censoring time from unif(0,200)
t=ifelse(t0<=ct, t0, ct)
delta=1*(t0<=ct)  ## about 75% event rate
ncase.train=sum(delta)
ncase.train / n_train ## event rate
summary(t) ## simulated survival time


#### non-linear 10 sine terms -----
method='nl-sine'
# set.seed(seed)
# err = matrix(rnorm(nrow(xtrue)*ncol(xtrue), 0, 0.25), ncol=nrow(xtrue)) 
interact.col = t(sin(t(xtrue.log)))
View(interact.col)
plot(xtrue.log["hsa.miR.20b..1",], interact.col["hsa.miR.20b..1",])

sd.interact = apply(interact.col, 1, sd)
interact.beta= rep(c(1,-1),5)*round(1/sd.interact/5, 2)
interact.beta

xtrue.nl3.log = interact.col
beta0.nl3 = interact.beta


h0=1e-2
t0=h=-log(runif(n_train)) * h0 * exp(t(xtrue.nl3.log)%*%beta0.nl3)

hist(t0)
set.seed(seed)
ct=runif(n_train, 0, 200) ## generate censoring time from unif(0,200)
t=ifelse(t0<=ct, t0, ct)
delta=1*(t0<=ct)  ## about 75% event rate
ncase.train=sum(delta)
print(paste0("Event rate: ", ncase.train / n_train))
summary(t)## simulated survival time 
# file.name = "simSurvival_nl-sine_10000_latest_orig.scale.csv"
file.name = "simSurvival_nl-sine_10000_latest_RW.scale.csv"

#### non-linear Gaussian hazard -----
method='nl-gauss'

n.selected = 10
set.seed(seed) # randomly select a interaction gene
selected.gene = sample(selected.geneid, n.selected)# "hsa.miR.887.1"
selected.gene

# scale gene expression to have stdev of 1
sel.xtrue.log.scaled = t(scale(t(xtrue.log[selected.gene,])))
x.mean= apply(sel.xtrue.log.scaled, 1, mean)
x.mean
# x.mean=read.csv('selected.gene.mean.scaled.csv')

# transform genes with a gaussian kernel 
gamma=1.5 # adjust this parameter for 
xtrue.nl4.scaled = exp(-gamma * ( (sel.xtrue.log.scaled)^2))
# View(interact.col.scaled)
# plot((sel.xtrue.log.scaled-x.mean)[1,], interact.col.scaled[1,])

sd.interact = apply(xtrue.nl4.scaled, 1, sd)
beta0.nl4= rep(c(1,-1),5)*round(1/sd.interact, 1)
beta0.nl4

h0 = 51
t0 = -log(runif(n_train)) * h0 / exp(t(xtrue.nl4.scaled)%*%beta0.nl4)
hist(t0)

ct=runif(n_train, 0, 200) ## generate censoring time from unif(0,200)
t=ifelse(t0<=ct, t0, ct)
delta=1*(t0<=ct)  ## about 75% event rate
ncase.train=sum(delta)
print(paste0("Event rate: ", ncase.train / n_train))
summary(t)## simulated survival time 

file.name ="simSurvival_nl-gauss_10000_101724_RW.scale.csv"

## Save survival data -----
out = data.frame(time=t, status=delta)
sim.survival = cbind(working, sim.survival.out)

out.dir  = file.path('data', dataType)
dir.create(out.dir, showWarnings = F)
write.csv(out, file=file.path(out.dir, file.name), row.names=F) ## save simulated survival outcomes
# write.csv(sim.survival,
#           file=file.path('data', dataType, paste0("simulate_survival_",method,"_",N*2,"_",date,".csv"))
# )


## K-M plot of simulated survival data -------
fit = survfit(Surv(time, status)~1, data=sim.survival, conf.type="log-log")
fig <- ggsurvfit(fit) +
  labs(
    title = 'Survival',
    x = 'Time (Days)',
    y = 'Survival Probability'
  ) +
  coord_cartesian(xlim = c(0, 200), ylim = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 200, by=20)) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2)) +
  add_confidence_interval() +
  add_risktable() +
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  theme_bw()
fig




# ------- START HERE ---------

# Function for quick simulation -----------
simulate_T <- function(h0, h.log, n_train, max.censor.time) {
  # simulate uncensored survival time
  t0=-log(runif(n_train)) *h0 / exp(h.log)
  
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



# Run Simulation -------------
# Linear exponential risk
linear.exp.rw.out = sim.survdata(surv.data=working, sim.dataType = 'linear', sim.method = 'exponential',
                                 h0=23, beta0=beta0.scaled , save_data = F, seed=1234,
                                 file.name="simSurvival_exponential_10000_091624_RW.scale.csv")
linear.exp.rw.surv= linear.exp.rw.out[[1]]
library(ggsurvfit)
library(survival)
plot(survfit(Surv(linear.exp.rw.surv$time,linear.exp.rw.surv$status)~1))

hist(linear.exp.rw.surv$time)

# 1 quadratic interaction
nl.1.rw.out = sim.survdata(surv.data = working, sim.dataType = 'nonlinear', sim.method = 'nl-1',
                           h0=750, beta0=beta0.scaled , save_data = T,
                           file.name="simSurvival_nl-1_10000_091624_RW.scale.csv")
nl.1.rw.surv=nl.1.rw.out[[1]]

# Quadratic interaction
nl.qua.rw.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method = 'nl-quadratic',
                             h0=2.5, beta0=beta0.scaled, save_data = T,
                             file.name="simSurvival_nl-quadratic_10000_091624_RW.scale.csv")
nl.qua.rw.surv = nl.qua.rw.out[[1]]

# Sine interaction
nl.sine.rw.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-sine',
                              h0=3.8, beta0=beta0.scaled, save_data = T,
                              file.name="simSurvival_nl-sine_10000_091624_RW.scale.csv")
nl.sine.rw.surv = nl.sine.rw.out[[1]]

# Gaussian kernel
nl.gauss.rw.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-gauss',
                               h0=51, simulate_beta = F,
                               save_data = F,
                               file.name="simSurvival_nl-gauss_10000_092624_RW.scale.csv")
nl.gauss.rw.surv = nl.gauss.rw.out[[1]]


# Run Simulation with unscaled betas -----
linear.exp.out = sim.survdata(surv.data=working, sim.dataType = 'linear', sim.method = 'exponential',
                              h0=0.4e-3, beta0=beta0, save_data = T, seed=1234,
                              file.name="simSurvival_exponential_10000_091824_orig.scale.csv")
linear.exp.surv= linear.exp.out[[1]]


nl.1.out = sim.survdata(surv.data = working, sim.dataType = 'nonlinear', sim.method = 'nl-1',
                        h0=1.7e-2, beta0=beta0, save_data = T, seed=1234,
                        file.name="simSurvival_nl-1_10000_091824_orig.scale.csv")
nl.1.surv=nl.1.out[[1]]


nl.qua.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method = 'nl-quadratic',
                          h0=4e-5, beta0=beta0, save_data = T,
                          file.name="simSurvival_nl-quadratic_10000_091824_orig.scale.csv")
nl.qua.surv = nl.qua.out[[1]]


nl.sine.out = sim.survdata(surv.data=working, sim.dataType='nonlinear', sim.method='nl-sine',
                           h0=1.5e-4, beta0=beta0, save_data = T,
                           file.name="simSurvival_nl-sine_10000_091824_orig.scale.csv")
nl.sine.surv = nl.sine.out[[1]]

View(cbind(linear.exp.surv, nl.1.surv, nl.qua.surv, nl.sine.surv))
