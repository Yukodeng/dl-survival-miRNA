#Run simulation with following command line:
#   sbatch --array=1-20 Batchsub<xxx>.sh 
#                 OR 
#   qsubR -q miniQ [-m mem] -a 1-20%10 [-s xxx.R]
#
#DATE       NOTES
#27JAN2025  Modified directory pathway 
#           Updated lines to add to the bash file to follow the configuration of the Biost omega00 cluster

############ Parametric simulation with mainSimSeqedgeRpar.R with batch effects on mean and observed counts (Nov 2024) #############
setwd("~/dl-survival-miRNA/simulation/batch_effect/Rscripts")
# setwd("C:/Users/ni.304/OneDrive - The Ohio State University/Research/R21 Nov2020 submission with Li-Xuan Qin/Augmented sequencing data analysis")

SubFileGen=function(beta0, nonzero_position, n_train, n_test, ps_rho_cutoff, beta_sort_train, beta_sort_test, batch_size, he_train, he_test, 
                    batch_sd_mu, batch_sd_obs, batch_sort_train, batch_sort_test, stratify_train, stratify_test, cmbt_train, cmbt_test, cmbt_frozen,
                    norm_type, norm_train, norm_test, standardize, univ_cutoff_range, nuniv_cutoff, lambda_max_glmnet, 
                    nlambda, nfold, addon, nsim, nodes, ppn, walltime){
  func.call=paste(#"source('/users/PAS1476/andyni/Augseq/mainSimSeqedgeRparV2.R')\n
                  "source('~/dl-survival-miRNA/simulation/batch_effect/mainSimSeqedgeRparV2.R')\n
seed <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))\n
mainSimSeqedgeRpar(", beta0, ",", nonzero_position, ",", n_train, ",", n_test, ",", ps_rho_cutoff, ",", 
                  beta_sort_train, ",", beta_sort_test, ",", batch_size, ",", he_train, ",", he_test, ",", batch_sd_mu, ",", batch_sd_obs,
                  ",", batch_sort_train, ",", batch_sort_test, ",", stratify_train, ",", stratify_test, ",", cmbt_train, ",", cmbt_test, ",", cmbt_frozen, ",", norm_type, ",", norm_train, ",", 
                  norm_test, ",", standardize, ",", univ_cutoff_range, ",", nuniv_cutoff, ",", lambda_max_glmnet, ",", nlambda, ",", nfold, ",", addon, ",", nsim, ", seed)", sep="")
  call.file.name=paste("Callhe", he_train*batch_sd_mu, he_train*batch_sd_obs, he_test*batch_sd_mu, he_test*batch_sd_obs, "norm", norm_type*norm_train, norm_type*norm_test, "tsort", beta_sort_train, beta_sort_test, ".R", sep="")
  write.table(func.call, call.file.name, row.names=F, col.names=F, quote=F)
  
  # cd /users/PAS1476/andyni/Augseq\n
  #SBATCH --account=PAS1476 \n
  batch.sub=paste("#!/bin/bash\n
#SBATCH --job-name=batch-sim\n
#SBATCH --output=batch-sim.log\n
#SBATCH --error=batch-sim_error.log\n
#SBATCH --partition=a2gpu\n  
#SBATCH --qos=a2gpu\n
#SBATCH --gres=gpu:A2:1\n
#SBATCH --nodes=", nodes, " --ntasks-per-node=", ppn, "\n
#SBATCH --time=", walltime, "\n
module load R/4.0.2-gnu9.1\n
cd simulation/batch_effect/Rscripts\n
cp ", call.file.name, " $TMPDIR\n
cd $TMPDIR\n
Rscript ", call.file.name, " $SLURM_ARRAY_TASK_ID ", call.file.name, "out \n
cp * $SLURM_SUBMIT_DIR", sep="")
  write.table(batch.sub, paste("Batchsubhe", he_train*batch_sd_mu, he_train*batch_sd_obs, he_test*batch_sd_mu, he_test*batch_sd_obs,
                               "norm", norm_type*norm_train, norm_type*norm_test, "tsort", beta_sort_train, beta_sort_test, ".sh", sep=""), row.names=F, col.names=F, quote=F)
}


# ######### Large gene effects ##########
# 
# ### no train HE, no test HE, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he0000bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")  ## need 20 jobs to make 400 replicates
# 
# ### no train HE, no test HE, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he0000bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # no train HE, no test HE, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he0000bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### no train HE, no test HE, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he0000bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### no train HE, no test HE, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he0000bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# 
# 
# 
# ### train HE, no test HE, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")  ## need 20 jobs to make 400 replicates
# 
# ### train HE, no test HE, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # train HE, no test HE, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, no test HE, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, no test HE, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# 
# 
# 
# ### train HE, no test HE, train sort, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm00tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")  ## need 20 jobs to make 400 replicates
# 
# ### train HE, no test HE, train sort, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm11tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # train HE, no test HE, train sort, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm22tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, no test HE, train sort, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm33tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, no test HE, train sort, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4400bs20norm44tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# 
# 
# 
# ### train HE, test HE, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")  ## need 20 jobs to make 400 replicates
# 
# ### train HE, test HE, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # train HE, test HE, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00") 
# 
# 
# 
# 
# ### train HE, test HE, train sort, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm00tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="8:00:00")  ## need 40 jobs to make 400 replicates
# 
# ### train HE, test HE, train sort, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm11tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # train HE, test HE, train sort, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm22tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, train sort, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm33tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, train sort, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm44tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# 
# 
# 
# ### train HE, test HE, test sort, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm00tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00")  ## need 40 jobs to make 400 replicates
# 
# ### train HE, test HE, test sort, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm11tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # train HE, test HE, test sort, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm22tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, test sort, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm33tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, test sort, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm44tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00")
# 
# 
# 
# 
# ### train HE, test HE, train sort, test sort, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm00tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00")  ## need 40 jobs to make 400 replicates
# 
# ### train HE, test HE, train sort, test sort, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm11tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # train HE, test HE, train sort, test sort, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm22tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, train sort, test sort, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm33tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, train sort, test sort, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm44tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00")
# 
# 
# 
# 
# ### train HE, test HE, train sort, test reverse sort, no norm, batch size 20, largest batch effect ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm00tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00")  ## need 40 jobs to make 400 replicates
# 
# ### train HE, test HE, train sort, test reverse sort, TC norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm11tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# # train HE, test HE, train sort, test reverse sort, upper quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm22tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, train sort, test reverse sort, TMM norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm33tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00") 
# 
# ### train HE, test HE, train sort, test reverse sort, quantile norm, batch size 20 ###
# SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
#            ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=4, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
#            norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
#            addon="'he4444bs20norm44tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="8:00:00")









###################################### Large gene effects, mean-level batch effects only ###################################

### no train HE, no test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="8:00:00")  ## need 20 jobs to make 400 replicates

### no train HE, no test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# no train HE, no test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### no train HE, no test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### no train HE, no test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, no test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, no test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, no test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, no test HE, train sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm00tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, no test HE, train sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm11tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, no test HE, train sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm22tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, train sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm33tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, train sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2000bs20norm44tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, train sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm00tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm11tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm22tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm33tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm44tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, test sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm00tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, test sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm11tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, test sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm22tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, test sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm33tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, test sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm44tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")




### train HE, test HE, train sort, test sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm00tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, test sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm11tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, test sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm22tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm33tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm44tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")




### train HE, test HE, train sort, test reverse sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm00tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, test reverse sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm11tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, test reverse sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm22tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test reverse sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm33tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test reverse sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2020bs20norm44tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")








######################################### Large gene effects, observation-level batch effects only ######################################

### no train HE, no test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### no train HE, no test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# no train HE, no test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### no train HE, no test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### no train HE, no test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, no test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, no test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, no test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, no test HE, train sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm00tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, no test HE, train sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm11tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, no test HE, train sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm22tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, train sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm33tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, train sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0200bs20norm44tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, train sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm00tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm11tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm22tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm33tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm44tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, test sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm00tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, test sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm11tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, test sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm22tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, test sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm33tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, test sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm44tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")




### train HE, test HE, train sort, test sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm00tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, test sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm11tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, test sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm22tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm33tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm44tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")




### train HE, test HE, train sort, test reverse sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm00tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, test reverse sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm11tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, test reverse sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm22tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test reverse sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm33tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test reverse sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=0, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0202bs20norm44tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")












######################################### Large gene effects, both mean-level and observation-level batch effects ######################################

### no train HE, no test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### no train HE, no test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# no train HE, no test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### no train HE, no test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### no train HE, no test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=0, he_test=0, batch_sd_mu=0, batch_sd_obs=0, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he0000bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, no test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm00tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, no test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm11tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, no test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm22tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm33tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm44tsort00'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, no test HE, train sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm00tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00")  ## need 20 jobs to make 400 replicates

### train HE, no test HE, train sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm11tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, no test HE, train sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm22tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, train sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm33tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, no test HE, train sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=0, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2200bs20norm44tsort0.050'", nsim=20, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm00tsort00'", nsim=5, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 200 replicates

### train HE, test HE, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm11tsort00'", nsim=5, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm22tsort00'", nsim=5, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm33tsort00'", nsim=5, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=0, cmbt_test=0, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm44tsort00'", nsim=5, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, train sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm00tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm11tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm22tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm33tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm44tsort0.050'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 




### train HE, test HE, test sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm00tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, test sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm11tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, test sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm22tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, test sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm33tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, test sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm44tsort00.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")




### train HE, test HE, train sort, test sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm00tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, test sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm11tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, test sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm22tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm33tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm44tsort0.050.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")




### train HE, test HE, train sort, test reverse sort, no norm, batch size 20, largest batch effect ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=0, norm_train=0, norm_test=0, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm00tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")  ## need 40 jobs to make 400 replicates

### train HE, test HE, train sort, test reverse sort, TC norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=1, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm11tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

# train HE, test HE, train sort, test reverse sort, upper quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=2, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm22tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test reverse sort, TMM norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=3, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm33tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00") 

### train HE, test HE, train sort, test reverse sort, quantile norm, batch size 20 ###
SubFileGen(beta0="c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)", nonzero_position="c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)", n_train=300, n_test=200,
           ps_rho_cutoff=0.9, beta_sort_train=0.05, beta_sort_test=-0.05, batch_size=20, he_train=1, he_test=1, batch_sd_mu=2, batch_sd_obs=2, batch_sort_train=5, batch_sort_test=5, stratify_train=1, stratify_test=1, cmbt_train=1, cmbt_test=1, cmbt_frozen=0, 
           norm_type=4, norm_train=1, norm_test=1, standardize=0, univ_cutoff_range="c(0,0.05)", nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=5,
           addon="'he2222bs20norm44tsort0.05-0.05'", nsim=10, nodes=1, ppn=1, walltime="16:00:00")

