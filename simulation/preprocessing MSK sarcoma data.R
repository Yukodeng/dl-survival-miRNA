############ use edgeR to estimate mu and phi for each gene using original MSK dataset (n=27+27)
library(edgeR)
library(dplyr)
library(reticulate)
# installing precision.seq
# install.packages("combinat")                     
# install.packages("https://cran.r-project.org/src/contrib/Archive/PoissonSeq/PoissonSeq_1.1.1.tar.gz")
# devtools::install_github("LXQin/precision.seq")

setwd("~/dl-survival-miRNA/")

# Load RWD ----------------------------------------------------------
# mirna.names = colnames(read.csv("MSKpair_MXF_merge.csv")[2:2067])
# mirna.names = as.data.frame(mirna.names)
# write.csv(mirna.names, "../data/MSKpair_miRNA_names.csv")

# mxf.rwd=read.csv("MSKpair_MXF_merge.csv")[,2:1034]
# pmfh.rwd=read.csv("MSKpair_PMFH_merge.csv")[,2:1034]
load(file.path('raw-data', 'precision.data.RData'))
workingt.rwd = rbind(data.benchmark, data.test)

workingt.rwd.log = log2(workingt.rwd+1)
rwd.max = max(workingt.rwd.log)
rwd.sd  = sd(workingt.rwd.log)
# rwd.max = max(workingt.rwd)
# rwd.sd  = sd(workingt.rwd)
print(rwd.max)
print(rwd.sd)

# ## Get maximums per gene ----
# working.rwd=rbind(mxf.rwd, pmfh.rwd)
# 
# # working.rwd %>% View()
# rwd.max=apply(working.rwd, 2, max)
# rwd.sd=apply(working.rwd, 2, sd)
# 
# plot(x = seq(1, length(rwd.sd[, gene.])), y = rwd.sd[rwd.sd<1e5])
# rwd.sd[rwd.sd>1e5]
# hist(working.rwd[,'hsa.miR.21.1'])
# write.csv(rwd.max, "~/dl-survival-miRNA/data/miRNA_max_count.csv")
# 
# surv.data[,"hsa.miR.21.1 "]

 
# Load augmented data with/without batch --------------------------------
### Attempt 11APR2025; N=10000/batch
N = 10000
P = 538 
n_batch = 10
g.names = read.csv("data/augmentation_miRNA_names_538.csv")$gene      # gene names

surv.data=data.frame(matrix(NA, nrow=0, ncol=P*2))
for(i in 1:n_batch){
  dat.sub=read.csv(
    file.path("raw-data", "with_without_batch_maf", paste0("batch-",i,"_epochES_maf_generated.csv")),
    header=F
  )[1:N,1:(P*2)]  ## test data n=200 from each batch
  
  cat(glue::glue("Batch {i}:"), sum(is.infinite(unlist(dat.sub))) / (N*P) * 100,"%\n") ## percentage of inf values in each batch
  # cat(glue::glue("Batch {i}:"), mean(apply(dat.sub, 1, function(row) any(is.infinite(row)))) * 100, "%\n") # percentage of samples w inf values in each batch
  
  ## merge batches
  colnames(dat.sub)=g.names
  surv.data = rbind(surv.data, dat.sub)
}


# Capping values
cap_and_add_noise <- function(x, cap, sd) {
  new_cap = cap + 0.1*sd
  xs = x[is.infinite(x) | (x > new_cap)]
  xs_hat = pmin( xs, new_cap + rnorm(length(xs), 0, 0.01*sd) )
  x[is.infinite(x) | (x > new_cap)] = xs_hat
  return(x)
}

surv.data.cap = (2 ^ as.data.frame(lapply(log2(surv.data+1), cap_and_add_noise, cap=rwd.max, sd=rwd.sd)) - 1) #%>% View()
max(surv.data.cap)


# Test reticulate -------------------------------------------------------
skl_ms <- import("sklearn.model_selection")

n = 500
n_run = 10

sim_train_py <- r_to_py(sim.train)
print(glue::glue("Running for N={n}..."))

for (run in 1:n_run) {
  run = 1
  if (n < nrow(sim.train)) {
    sim_train_py$strata <- sim_train_py$status$astype("str")$str$cat(
      sim_train_py[['batch.id']]$astype("str"), sep = "_"
    )
    
    train_sub_py = skl_ms$train_test_split(
      sim_train_py,
      train_size=as.integer(n), # train_size= n/self.train_df.shape[0], 
      shuffle=TRUE,
      random_state=as.integer(run),
      stratify = sim_train_py$strata
    )[[1]]
  } else {
    train_sub_py = sim_train_py
  }
  
  data_train = py_to_r(train_sub_py)
  data_test = sim.test
  x0_train = data_train[,selected.geneid]
  x0_test  = data_test[,selected.geneid]
  
  # oracle analysis -----------
  coxfit_o = tryCatch(
      coxph(Surv(data_train$time, data_train$status) ~ as.matrix(x0_train)), 
      error=function(e) e,
      warning=function(w) w
      )
  if(is(coxfit_o, "warning") | is(coxfit_o, "error")) {
    coxfit_o = coxph(Surv(data_train$time, data_train$status) ~ ridge(as.matrix(x0_train), theta=1))
  }
  c_o[run]=summary(coxfit_o)$concordance[1]      
  
  coxfit_o_test = coxph(
    Surv(data_test$time, data_test$status) ~ as.matrix(x0_test), 
    init=summary(coxfit_o)$coefficient[,1], 
    control=coxph.control(iter.max=0)
    )
  c_o_test[run]=summary(coxfit_o_test)$concordance[1]   
  
}
