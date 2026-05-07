############################################################################# ##
# Date      Note
# 15APR2026 Change to run all linear/nl marker sizes per ML model in 1 bash file  
############################################################################# ## 
library(dplyr)
library(glue)

# ################################## Function ##################################

generate_subFile <- function(sim.dataType = NULL, # E.g., "linear", "nl-quadratic"
                          p_vec = c(2,5,10,30),
                          he_train = 0,
                          he_test = 0,
                          beta_sort_train = 0,
                          beta_sort_test = 0,
                          norm_type = 0,
                          keywords = NULL,
                          subset_sizes = c(100,200,500,1000,2000,5000,10000),
                          runs_per_size = c(20,20,20,20,20,20,20),
                          splits_per_size = c(3,5,5,5,10,10,10),
                          trials_per_size = c(30,30,30,50,50,50,50),
                          trial_threshold = 30,
                          GPUtype = 'l4') {

  norm_map <- c('None','TC','UQ','TMM','Quantile','DEseq','Med')
  dataTypes_vec <- paste0(sim.dataType,"-p", p_vec)
  
  dataTypes_str <- paste(dataTypes_vec, collapse = " ")
  subsets_str <- paste(subset_sizes, collapse = ",")
  runs_str   <- paste(runs_per_size, collapse = ",")
  splits_str <- paste(splits_per_size, collapse = ",")

  # Code scenario name -----------------------------------------------
  
  ## E.g., "BE11Asso00_normTC" ({batchType}_norm{normType})
  convert_num_to_indicator <- function(x) {
    case_when(
      x>0  ~ 1,
      x==0 ~ 0,
      x<0  ~ -1
    )
  }
  sort_train = convert_num_to_indicator(beta_sort_train)
  sort_test =  convert_num_to_indicator(beta_sort_test)
  batchNormType = glue(
    "BE{he_train}{he_test}Asso{sort_train}{sort_test}_norm{norm_map[norm_type+1]}"
  )

  # Generate Slurm script for job submission -------------------------
  
  for (dataType in dataTypes_vec) {

    ## DL config file --------------------------------------
    dir.create(here::here('configs', batchNormType), recursive=T, showWarnings = F)

    config = list(
      batchNormType = batchNormType,
      dataName = dataType,
      storage_url = "sqlite:///deepsurv-torch-hp-log.db",
      keywords = keywords,
      time_col = "time",
      status_col = "status",
      batch_col = "batch_id",
      keep_batch = TRUE,
      hyperparameters = list(
        num_nodes = list(
          type = "categorical",
          choices = list(list(128), list(64), list(32), #list(16),
                        c(128,64), c(64,32), c(32,16), 
                        c(64,64,32), c(32,32,16))),
        dropout = list(type = "float", low = 0.1, high = 0.5),
        learning_rate = list(type = "float", low = 1e-4, high = 5e-3, log = TRUE),
        weight_decay = list(type = "float", low = 1e-4, high = 5e-2, log = TRUE),
        batch_size = list(type = "categorical", choices = c(32, 64, 128))
      ),
      subset_sizes = subset_sizes,
      runs_per_size = runs_per_size,
      splits_per_size = splits_per_size,
      trials_per_size = trials_per_size,
      trial_threshold = trial_threshold,
      n_jobs = 1,
      is_tune = TRUE,
      is_save = TRUE,
      early_stop_per_size = list(
        `100` = list(patience = 20, min_delta = 5e-3),
        `200` = list(patience = 20, min_delta = 5e-3),
        `500` = list(patience = 20, min_delta = 5e-3),
        `1000` = list(patience = 30, min_delta = 1e-3),
        `2000` = list(patience = 30, min_delta = 5e-4),
        `5000` = list(patience = 35, min_delta = 1e-4),
        `10000` = list(patience = 40, min_delta = 1e-4)
      )
    )
    # Save to json files
    jsonlite::write_json(config, 
      here::here('configs', batchNormType, glue("{batchNormType}-{dataType}.json")), pretty=T, auto_unbox=T)
    jsonlite::write_json( # [Added 08/11/25] Added stratified DeepSurv configuration file
      append(config, list(stratified = TRUE)),
      here::here('configs', batchNormType, glue("{batchNormType}-{dataType}-stratified.json")),
      pretty=T, auto_unbox=T
    )
  
    ## DL Slurm script --------------------------------------

    ## [UPDATE 09052025] updated gres=gpu: to gres=shard: per Venkat suggestion
    bash.sub = glue("#!/bin/bash
#SBATCH --job-name={batchNormType}-{dataType}{GPUtype}
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition={tolower(GPUtype)}gpu
#SBATCH --qos={tolower(GPUtype)}gpu
#SBATCH --gres=shard:{toupper(GPUtype)}:1
#SBATCH --time=12:00:00
exec > >(tee -a jobs/deepsurv/logs/{batchNormType}/{batchNormType}-{dataType}{GPUtype}.log) 2>&1

export CUDA_HOME=/opt/cuda118
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$CUDA_HOMSE/bin:$PATH\n
source /home/nfs/dengy/dl-surv/bin/activate
pss-deepsurv --config configs/{batchNormType}/{batchNormType}-{dataType}.json
pss-deepsurv --config configs/{batchNormType}/{batchNormType}-{dataType}-stratified.json"
    )
    
    dir.create(here::here('jobs', 'deepsurv', batchNormType), recursive=T, showWarnings = F)
    dir.create(here::here('jobs', 'deepsurv', 'logs', batchNormType), recursive=T, showWarnings = F)
    write.table(bash.sub, 
      here::here('jobs','deepsurv', batchNormType, glue("{batchNormType}-{dataType}.sh")), 
      row.names=F, col.names=F, quote=F)
  }
  
  ## ML bash script ---------------------------------------

  bash.sksurv.sgb <- glue(
"#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}-sgb
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G

LOGFILE='jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-sgb.log'
exec > >(tee -a $LOGFILE) 2>&1

source ~/dl-surv/bin/activate

COMMON_ARGS='--batchNormType {batchNormType} --modelType sgb --is_save \\
--subsets {subsets_str} \\
--runs {runs_str} \\
--splits {splits_str}'

for DATA_NAME in {dataTypes_str}; do
  echo \"==================================================\"
  echo \"Running SGB for $DATA_NAME (non-stratified)\"
  pss-sksurv --dataName $DATA_NAME $COMMON_ARGS
  echo \"Running SGB for $DATA_NAME (stratified)\"
  pss-sksurv --dataName $DATA_NAME $COMMON_ARGS --is_stratified
done

echo 'Congrats! All survival gradient boosting (SGB) model runs completed.'")

  ### SSVM ---------------
  
  bash.sksurv.ssvm <- glue(
"#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}-ssvm
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G

LOGFILE='jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-ssvm.log'
exec > >(tee -a $LOGFILE) 2>&1

source ~/dl-surv/bin/activate

COMMON_ARGS='--batchNormType {batchNormType} --modelType ssvm --is_save \\
--subsets {subsets_str} \\
--runs {runs_str} \\
--splits {splits_str}'

for DATA_NAME in {dataTypes_str}; do
  echo \"==================================================\"
  echo \"Running SSVM for $DATA_NAME (non-stratified)\"
  pss-sksurv --dataName $DATA_NAME $COMMON_ARGS
  echo \"Running SSVM for $DATA_NAME (stratified)\"
  pss-sksurv --dataName $DATA_NAME $COMMON_ARGS --is_stratified
done

echo 'Congrats! All survival support vector machine (SSVM) model runs completed.'")
  
  ### RSF ---------------
  
  bash.sksurv.rsf <- glue("#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}-rsf
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=64G
LOGFILE='jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-rsf.log'
exec > >(tee -a $LOGFILE) 2>&1

source ~/dl-surv/bin/activate

COMMON_ARGS='--batchNormType {batchNormType} --is_save --modelType rsf \\
--subsets 100,200,500,1000,2000,5000 \\
--runs 20,20,20,20,20,20 \\
--splits 3,5,5,5,10,10'

for DATA_NAME in {dataTypes_str}; do
  echo \"==================================================\"
  echo \"Running RSF for $DATA_NAME (non-stratified)\"
  pss-sksurv --dataName $DATA_NAME $COMMON_ARGS
  echo \"Running RSF for $DATA_NAME (stratified)\"
  pss-sksurv --dataName $DATA_NAME $COMMON_ARGS --is_stratified
done

echo 'Congrats! All random survival forest (RSF) model runs completed.'")
  
  ### RSF (8000) ---------------------------------
  
  bash.sksurv.rsf.8000 <- glue("#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}-8000
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --array=0-4%1
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
exec > >(tee -a jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-rsf-8000.log) 2>&1

source ~/dl-surv/bin/activate

IDX=$SLURM_ARRAY_TASK_ID
SEEDS=(0 1 2 3 4)
FILE=\"model_results_8000.csv\"

pss-sksurv --batchNormType {batchNormType} --dataName {sim.dataType} \\
--subsets 8000 --runs 1 --seeds ${{SEEDS[$IDX]}} --fileName ${{FILE}} \\
--is_save --modelType rsf
pss-sksurv --batchNormType {batchNormType} --dataName {sim.dataType} \\
--subsets 8000 --runs 1 --seeds ${{SEEDS[$IDX]}} --fileName ${{FILE}} \\
--is_save --modelType rsf --is_stratified

echo 'Finished all RSF (N=8000) model runs.'")
  
  # [9/18/2025] save slurm script for GB models separately
  file.ssvm = here::here('jobs', 'sksurv', batchNormType, glue("{batchNormType}-{sim.dataType}-ssvm.sh"))
  file.sgb  = here::here('jobs', 'sksurv', batchNormType, glue("{batchNormType}-{sim.dataType}-sgb.sh"))
  file.rsf  = here::here('jobs', 'sksurv', batchNormType, glue("{batchNormType}-{sim.dataType}-rsf.sh"))
  # file.rsf.8000 = here::here('jobs', 'sksurv', batchNormType, glue("{batchNormType}-{sim.dataType}-rsf-8000.sh"))
  
  dir.create(here::here('jobs', 'sksurv', batchNormType), recursive = TRUE, showWarnings = FALSE)
  dir.create(here::here('jobs', 'sksurv', 'logs', batchNormType), recursive = TRUE, showWarnings = FALSE)

  writeLines(bash.sksurv.sgb, file.sgb)
  writeLines(bash.sksurv.ssvm, file.ssvm)
  writeLines(bash.sksurv.rsf, file.rsf)
  # writeLines(bash.sksurv.rsf.8000, file.rsf.8000)

  cat('Done.\n+++++++++++++++++++++++++++++++++++++++++++++++\n')
}



# ############################## Bash Parameters ##################################

sim.dataType.vec <- c(
  "linear",
  "nl-quadratic",
  "nl-interaction",
  "nl-sine"
)
p_vec <- c(2, 5, 10, 30)
# dataType_combs <- expand.grid(dataType_vec, p_vec, stringsAsFactors = FALSE)
# dataType_list <- paste(dataType_combs$Var1, dataType_combs$Var2, sep = "-p")
subset_sizes = c(100,200,500,1000,2000,5000)
runs_per_size = c(20,20,20,20,20,20)
splits_per_size = c(3,5,5,5,10,10)
trials_per_size = c(30,30,30,50,50,50)

## ############################# BE00Asso00 ###############################

date = NULL#"061825"

for (norm_type in 0:6) {

  purrr::map(sim.dataType.vec,
    ~ generate_subFile(
      sim.dataType = .x,
      p_vec = p_vec,
      he_train = 0,
      he_test = 0,
      beta_sort_train = 0,
      beta_sort_test = 0,
      norm_type = norm_type,
      subset_sizes = subset_sizes,
      runs_per_size = runs_per_size,
      splits_per_size = splits_per_size,
      trials_per_size = trials_per_size,
      GPUtype = 'l4')
  )
}
  

## ############################# BE10Asso00 ###############################

date = NULL #"071625"

for (norm_type in 0:6) {

  purrr::map(sim.dataType.vec,
    ~ generate_subFile(
      sim.dataType = .x,
      p_vec = p_vec,
      he_train = 1,
      he_test = 0,
      beta_sort_train = 0,
      beta_sort_test = 0,
      norm_type = norm_type,
      subset_sizes = subset_sizes,
      runs_per_size = runs_per_size,
      splits_per_size = splits_per_size,
      trials_per_size = trials_per_size,
      GPUtype = 'l4')
  )
}


## ############################# BE11Asso00 ###############################

keywords = NULL#"080425"

for (norm_type in 0:6) {

  purrr::map(sim.dataType.vec,
    ~ generate_subFile(
      sim.dataType = .x,
      p_vec = p_vec,
      he_train = 1,
      he_test = 1,
      beta_sort_train = 0,
      beta_sort_test = 0,
      norm_type = norm_type,
      subset_sizes = subset_sizes,
      runs_per_size = runs_per_size,
      splits_per_size = splits_per_size,
      trials_per_size = trials_per_size,
      GPUtype = 'l4')
  )
}

# ## ############################# BE10Asso10 ###############################

# date = "072125"

# for (norm_type in 0:6) {

#   # Linear risk with moderate effects
#   generate_subFile(sim.dataType='linear-moderate',
#                    he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Linear risk with weak effects
#   generate_subFile(sim.dataType='linear-weak',
#                    he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Squared terms
#   generate_subFile(sim.dataType='nl-quadratic',
#                    he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Squared terms with intercept
#   generate_subFile(sim.dataType='nl-shiftquad',
#                    he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Gene-gene interaction terms
#   generate_subFile(sim.dataType='nl-interaction',
#                    he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Sine interaction
#   generate_subFile(sim.dataType='nl-sine',
#                    he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')
# }


# ## ############################# BE11Asso10 ###############################
# date = "091625"

# for (norm_type in 0:6) {
  
#   # Linear risk with moderate effects
#   generate_subFile(sim.dataType='linear-moderate',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')
  
#   # Linear risk with weak effects
#   generate_subFile(sim.dataType='linear-weak',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')
  
#   # Squared terms
#   generate_subFile(sim.dataType='nl-quadratic',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')
  
#   # Squared terms with intercept
#   generate_subFile(sim.dataType='nl-shiftquad',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')
  
#   # Gene-gene interaction terms
#   generate_subFile(sim.dataType='nl-interaction',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')
  
#   # Sine interaction
#   generate_subFile(sim.dataType='nl-sine',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                    date=date, GPUtype='l4')
# }




# ## ############################# BE11Asso11 ###############################
# date = "091725"

# for (norm_type in 0:6) {

#   # Linear risk with moderate effects
#   generate_subFile(sim.dataType='linear-moderate',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Linear risk with weak effects
#   generate_subFile(sim.dataType='linear-weak',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Squared terms
#   generate_subFile(sim.dataType='nl-quadratic',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Squared terms with intercept
#   generate_subFile(sim.dataType='nl-shiftquad',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Gene-gene interaction terms
#   generate_subFile(sim.dataType='nl-interaction',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                    date=date, GPUtype='l4')

#   # Sine interaction
#   generate_subFile(sim.dataType='nl-sine',
#                    he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                    date=date, GPUtype='l4')
# }

