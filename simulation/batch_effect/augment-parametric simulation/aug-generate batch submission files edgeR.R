############ Application of parametric analysis on augmented batch data with aug-mainSimSeqedgeRpar.R  ######
#Run simulation with following command line:
#   sbatch --array=1-20 Batchsub<xxx>.sh 
#                 OR 
#   qsubR -q miniQ [-m mem] -a 1-20%10 [-s xxx.R]
#
#DATE       NOTES
#17Feb2025  Copied from job submission generation file for parametric method


SubFileGen=function(beta0, nonzero_position, n_train, n_test, beta_sort_train, beta_sort_test, 
                    train_batch_size, test_batch_size, n_batch, he_train, he_test, batch_sd_mu, batch_sd_obs, batch_sort_train, batch_sort_test, 
                    stratify_train, stratify_test, cmbt_train, cmbt_test, cmbt_frozen,
                    norm_type, norm_train, norm_test, standardize, ps_rho_cutoff, 
                    univ_cutoff_range, nuniv_cutoff, lambda_max_glmnet, nlambda, nfold,
                    addon, nsim, nodes, ppn, walltime){
  setwd("~/dl-survival-miRNA/simulation/batch_effect/Rscripts")
  
  func.call=paste0("source('~/dl-survival-miRNA/simulation/batch_effect/augment-parametric simulation/aug-mainSimSeqedgeRparV2.R')\n
seed <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))\n
mainSimSeqedgeRpar.aug(", beta0, ",", nonzero_position, ",", n_train, ",", n_test, ",",
                  beta_sort_train, ",", beta_sort_test, ",", train_batch_size, ",", test_batch_size, ",", n_batch, ",", 
                  he_train, ",", he_test, ",", batch_sd_mu, ",", batch_sd_obs, ",", batch_sort_train, ",", batch_sort_test, ",", stratify_train, ",", stratify_test, ",", cmbt_train, ",", cmbt_test, ",", cmbt_frozen, ",",
                  norm_type, ",", norm_train, ",", norm_test, ",", standardize, ",", ps_rho_cutoff, ",", univ_cutoff_range, ",", nuniv_cutoff, ",", lambda_max_glmnet, ",", nlambda, ",", nfold, ",", 
                  addon, ",", nsim, ",seed)")
  call.file.name=paste0("Callhe", he_train*batch_sd_mu, he_train*batch_sd_obs, he_test*batch_sd_mu, he_test*batch_sd_obs, "norm", norm_type*norm_train, norm_type*norm_test, "tsort", beta_sort_train, beta_sort_test, "-aug.R")
  print(call.file.name)
  write.table(func.call, call.file.name, row.names=F, col.names=F, quote=F)
 
#   batch.sub=paste("#!/bin/bash\n
# #SBATCH --job-name=batch-sim\n
# #SBATCH --output=batch-sim.log\n
# #SBATCH --error=batch-sim_error.log\n
# #SBATCH --partition=a2gpu\n  
# #SBATCH --qos=a2gpu\n
# #SBATCH --gres=gpu:A2:1\n
# #SBATCH --nodes=", nodes, " --ntasks-per-node=", ppn, "\n
# #SBATCH --time=", walltime, "\n
# module load R/4.0.2-gnu9.1\n
# cd simulation/batch_effect/Rscripts\n
# cp ", call.file.name, " $TMPDIR\n
# cd $TMPDIR\n
# Rscript ", call.file.name, " $SLURM_ARRAY_TASK_ID ", call.file.name, "out \n
# cp * $SLURM_SUBMIT_DIR", sep="")
  # write.table(batch.sub, paste("Batchsubhe", he_train*batch_sd_mu, he_train*batch_sd_obs, he_test*batch_sd_mu, he_test*batch_sd_obs,
  #                              "norm", norm_type*norm_train, norm_type*norm_test, "tsort", beta_sort_train, beta_sort_test, ".sh", sep=""), row.names=F, col.names=F, quote=F)
}
c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)

###################################### Large gene effects, mean-level batch effects only ###################################
########### no train HE and no test HE #############
### no train HE, no test HE, no norm, batch size 10, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=0, he_test=0, 
           batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, ps_rho_cutoff=0.9, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm00tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")  ## need 20 jobs to make 400 replicates

### no train HE, no test HE, TC norm, batch size 10 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=0, he_test=0, 
           batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5,stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, ps_rho_cutoff=0.9, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm11tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 

# no train HE, no test HE, upper quantile norm, batch size 10 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=0, he_test=0, 
           batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, ps_rho_cutoff=0.9, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm22tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 

### no train HE, no test HE, TMM norm, batch size 10 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=0, he_test=0, 
           batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, ps_rho_cutoff=0.9, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm33tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 

### no train HE, no test HE, quantile norm, batch size 10 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=0, he_test=0, 
           batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, ps_rho_cutoff=0.9, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm44tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 



########### Train HE and no test HE ###################
### train HE, no test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=1, he_test=0,
           batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0,
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm00tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")  ## need 20 jobs to make 400 replicates

### train HE, no test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=1, he_test=0, 
           batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0,
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm11tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")

# train HE, no test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=1, he_test=0, 
           batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0,
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm22tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")

### train HE, no test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=1, he_test=0, 
           batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0,
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm33tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")

### train HE, no test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, train_batch_size=30, test_batch_size=20, n_batch=10, he_train=1, he_test=0,
           batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0,
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm44tsort00-aug'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")





################# Directly run jobs ###################
source('~/dl-survival-miRNA/simulation/batch_effect/augment-parametric simulation/aug-mainSimSeqedgeRparV2.R')

##################### Large gene effects, mean-level batch effects only ###################################
### no train HE, no test HE, no norm, batch size 20, largest batch effect ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,0,0,0,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm00tsort00-aug',20,seed)

### no train HE, no test HE, TC norm, batch size 20, largest batch effect ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,1,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm11tsort00-aug',20,seed)

# no train HE, no test HE, upper quantile norm, batch size 10 ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,2,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm22tsort00-aug',20,seed)

### no train HE, no test HE, TMM norm, batch size 10 ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,3,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm33tsort00-aug',20,seed)

### no train HE, no test HE, quantile norm, batch size 10 ###
### NOTE: This scenario always produces error: "Error. Skipping." for all iterations
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,4,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm44tsort00-aug',20,seed)


############################# train HE, no test HE #############################
### train HE, no test HE, no norm ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,0,0,0,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm00tsort00-aug',20,seed)

### train HE, no test HE, TC norm, batch size 20 ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,1,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm11tsort00-aug',20,seed)

### train HE, no test HE, upper quantile norm, batch size 20 ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,2,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm22tsort00-aug',20,seed)

### train HE, no test HE, TMM norm, batch size 20 ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,3,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm33tsort00-aug',20,seed)

### train HE, no test HE, quantile norm, batch size 20 ###
mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518),300,200,0,0,30,20,10,1,0,2,0,5,5,1,1,0,0,0,4,1,1,0,0.9,c(0,0.05),20,0.3,30,5,'he2000bs20norm44tsort00-aug',20,seed)

