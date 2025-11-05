library(dplyr)
library(glue)

# ################################## Function ##################################

generate_subFile <- function(sim.dataType=NULL, # E.g., "linear-moderate", "nl-quadratic"
                         he_train=0,
                         he_test=0,
                         beta_sort_train=0,
                         beta_sort_test=0,
                         norm_type=0,
                         date = format(Sys.Date(), "%m%d%y"),
                         subset_sizes = c(100,500,1000,2000,5000,10000),
                         runs_per_size = c(20,20,20,20,20,20),
                         splits_per_size = c(3,5,5,10,10,10), 
                         GPUtype='l4') {
  
  
  norm_map = c('None','TC','UQ','TMM','Quantile','DEseq','Med')
  
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
  
  ## DL config file --------------------------------------
  
  config = list(
    batchNormType = batchNormType,
    dataName = sim.dataType,
    storage_url = "sqlite:///deepsurv-torch-hp-log-1.db",
    keywords = list(glue("{date}")),
    time_col = "time",
    status_col = "status",
    batch_col = "batch.id",
    keep_batch = TRUE,
    hyperparameters = list(
      num_nodes = list(
        type = "categorical",
        choices = list(list(128), list(64), list(32), #list(16),
                       c(128,64), c(64,32), c(32,16), 
                       c(64,64,32),c(32,32,16))),
      dropout = list(type = "float", low = 0.1, high = 0.5),
      learning_rate = list(type = "float", low = 1e-4, high = 5e-3, log = TRUE),
      weight_decay = list(type = "float", low = 1e-4, high = 5e-2, log = TRUE),
      batch_size = list(type = "categorical", choices = c(32, 64, 128))
    ),
    subset_sizes = c(100, 500, 1000, 2000, 5000, 10000),
    runs_per_size = c(20, 20, 20, 20, 20, 20),
    splits_per_size = c(3, 5, 5, 10, 10, 10),
    trials_per_size = c(30, 30, 30, 30, 30, 30),
    trial_threshold = 20,
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
  
  # [UNCOMMENT!] ----
  # # Save to json file
  # dir.create(here::here('configs', batchNormType), showWarnings = F)
  # 
  # file = here::here('configs', batchNormType, glue("{batchNormType}-{sim.dataType}.json"))
  # jsonlite::write_json(config, file, pretty=T, auto_unbox=T)
  # 
  # # [Added 08/11/25] Added stratified DeepSurv configuration file
  # jsonlite::write_json(append(config, list(stratified = TRUE)),
  #                      here::here('configs', batchNormType, glue("{batchNormType}-{sim.dataType}-stratified.json")),
  #                      pretty=T, auto_unbox=T)
  
  
  ## DL Slurm script --------------------------------------
  ## [UPDATE 09052025] updated gres=gpu: to gres=shard: per Venkat suggestion
  bash.sub = glue(
    "#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}{GPUtype}
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition={tolower(GPUtype)}gpu
#SBATCH --qos={tolower(GPUtype)}gpu
#SBATCH --gres=shard:{toupper(GPUtype)}:1
#SBATCH --time=12:00:00
exec > >(tee -a jobs/deepsurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}{GPUtype}.log) 2>&1

export CUDA_HOME=/opt/cuda118
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$CUDA_HOMSE/bin:$PATH\n
source /home/nfs/dengy/dl-surv/bin/activate
python3 main.py --config configs/{batchNormType}/{batchNormType}-{sim.dataType}.json
python3 main.py --config configs/{batchNormType}/{batchNormType}-{sim.dataType}-stratified.json\n
deactivate"
  )
  
  # [UNCOMMENT!] ----
  dir.create(here::here('jobs', 'deepsurv', batchNormType), showWarnings = F)
  dir.create(here::here('jobs', 'deepsurv', 'logs', batchNormType), showWarnings = F)
  
  # file = here::here('jobs','deepsurv', batchNormType, glue("{batchNormType}-{sim.dataType}.sh"))
  # write.table(bash.sub, file, row.names=F, col.names=F, quote=F)

  
  ## ML bash script ---------------------------------------
  
  ### SVM / RSF(N<=2000) ---------------
  
  bash.sksurv = glue("#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=12G
LOGFILE='jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}.log'
exec > >(tee -a $LOGFILE) 2>&1

source ~/dl-surv/bin/activate

SVM_COMMON_ARGS='--batchNormType {batchNormType} --dataName {sim.dataType} --keywords {date} --is_save \\
--subsets 100,500,1000,2000,5000,10000 \\
--runs 20,20,20,20,20,20 \\
--splits 3,5,5,10,10,10'

RSF_COMMON_ARGS='--batchNormType {batchNormType} --dataName {sim.dataType} --keywords {date} --is_save \\
--subsets 100,500,1000,2000 \\
--runs 20,20,20,20 \\
--splits 3,5,5,10'

echo \"==========================================\"
echo \"Running SVM (stratified) Models...\"
echo \"==========================================\"
# python main-sksurv.py $SVM_COMMON_ARGS --modelType svm
python main-sksurv.py $SVM_COMMON_ARGS --is_stratified --modelType svm

echo \"==========================================\"
echo \"Running RSF (stratified) Models...\"
echo \"==========================================\"
#python main-sksurv.py $RSF_COMMON_ARGS --modelType rsf
python main-sksurv.py $RSF_COMMON_ARGS --is_stratified --modelType rsf

echo 'Congrats! All survival models completed.'")
  
  
  ### GB ---------------
  
  bash.sksurv.gb = glue("#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}-gb
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=12G
LOGFILE='jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-gb.log'
exec > >(tee -a $LOGFILE) 2>&1

source ~/dl-surv/bin/activate

COMMON_ARGS='--batchNormType {batchNormType} --dataName {sim.dataType} --keywords {date} --is_save \\
--subsets 100,500,1000,2000,5000,10000 \\
--runs 20,20,20,20,20,20 \\
--splits 3,5,5,10,10,10'

# python main-sksurv.py $COMMON_ARGS --modelType gb
python main-sksurv.py $COMMON_ARGS --is_stratified --modelType gb

echo 'Congrats! All survival GB models completed.'")
  
  
  ### RSF (5000) -------------------------------
  
  slurm_arrayID = "0-4%1"
  slurm_string = "RUNS=(2 2 2 2 2)
SEEDS=(\"0,1\" \"2,3\" \"4,5\" \"6,7\" \"8,9\")"
  
  bash.sksurv.5000 <- glue("#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}-5000
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --array={slurm_arrayID}
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=36G
exec > >(tee -a jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-5000.log) 2>&1

source ~/dl-surv/bin/activate

IDX=$SLURM_ARRAY_TASK_ID
{slurm_string}
FILE=\"model_results_5000.csv\"


# python main-sksurv.py --batchNormType {batchNormType} --dataName {sim.dataType} \\
--subsets 5000 --runs ${{RUNS[$IDX]}} --seeds ${{SEEDS[$IDX]}} \\
--fileName ${{FILE}} \\
--keywords {date} \\
--is_save --modelType rsf
python main-sksurv.py --batchNormType {batchNormType} --dataName {sim.dataType} \\
--subsets 5000 --runs ${{RUNS[$IDX]}} --seeds ${{SEEDS[$IDX]}} \\
--fileName ${{FILE}} \\
--keywords {date} \\
--is_save --is_stratified --modelType rsf

echo 'Finished running all seeds.")
  
  
  ### RSF (8000) ---------------------------------
  
#   bash.sksurv.8000 <- glue("#!/bin/bash
# LOGFILE='jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-8000.log'
# exec > >(tee -a $LOGFILE) 2>&1
# 
# source ~/dl-surv/bin/activate
# 
# for SEED in {{0..4}}; do
#     echo \"==========================================\"
#     echo \"Running script with seed: $SEED...\"
#     echo \"==========================================\"
# 
#     python main-sksurv.py \\
# --batchNormType {batchNormType} \\
# --dataName {sim.dataType} \\
# --subsets 8000 \\
# --runs 1 \\
# --seeds $SEED \\
# --fileName \"model_results_8000.csv\" \\
# --keywords {date} \\
# --is_save \\
# --modelType rsf
# done
# 
# echo 'Finished running all seeds.'")
  
    bash.sksurv.8000 <- glue("#!/bin/bash
#SBATCH --job-name={batchNormType}-{sim.dataType}-8000
#SBATCH --error=slurm-temp.log
#SBATCH --output=slurm-temp.log
#SBATCH --array=0-4%1
#SBATCH --partition=epycQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=72G
exec > >(tee -a jobs/sksurv/logs/{batchNormType}/{batchNormType}-{sim.dataType}-8000.log) 2>&1

source ~/dl-surv/bin/activate

IDX=$SLURM_ARRAY_TASK_ID
SEEDS=(0 1 2 3 4)
FILE=\"model_results_8000.csv\"


# python main-sksurv.py --batchNormType {batchNormType} --dataName {sim.dataType} \\
--subsets 8000 --runs 1 --seeds ${{SEEDS[$IDX]}} \\
--fileName ${{FILE}} \\
--keywords {date} \\
--is_save --modelType rsf
python main-sksurv.py --batchNormType {batchNormType} --dataName {sim.dataType} \\
--subsets 8000 --runs 1 --seeds ${{SEEDS[$IDX]}} \\
--fileName ${{FILE}} \\
--keywords {date} \\
--is_save --is_stratified --modelType rsf

echo 'Finished running all seeds.'")
  
  
  # [9/18/2025] save slurm script for GB models separately
  file = here::here('jobs','sksurv',batchNormType,glue("{batchNormType}-{sim.dataType}.sh"))
  file.gb = here::here('jobs','sksurv',batchNormType,glue("{batchNormType}-{sim.dataType}-gb.sh"))
  file.5000 = here::here('jobs','sksurv',batchNormType,glue("{batchNormType}-{sim.dataType}-5000.sh"))
  file.8000 = here::here('jobs','sksurv',batchNormType,glue("{batchNormType}-{sim.dataType}-8000.sh")) # file.10000 = here::here('jobs','sksurv',glue("{batchNormType}-{sim.dataType}-10000.sh"))
  
  dir.create(here::here('jobs', 'sksurv', batchNormType), showWarnings = F)
  dir.create(here::here('jobs', 'sksurv', 'logs', batchNormType), showWarnings = F)
  
  # [UNCOMMENT!] ----
  write.table(bash.sksurv, file, row.names=F, col.names=F, quote=F)
  write.table(bash.sksurv.gb, file.gb, row.names=F, col.names=F, quote=F)
  write.table(bash.sksurv.5000, file.5000, row.names=F, col.names=F, quote=F)
  write.table(bash.sksurv.8000, file.8000, row.names=F, col.names=F, quote=F) # write.table(bash.sksurv.10000, file.10000, row.names=F, col.names=F, quote=F)


  cat('Done.\n+++++++++++++++++++++++++++++++++++++++++++++++\n')
}



# ############################## Bash Generation ##################################


## ############################# BE00Asso00 ###############################

date = "061825"

for (norm_type in 0:6) {
  
  # Linear risk with moderate effects
  generate_subFile(sim.dataType='linear-moderate',
                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Linear risk with weak effects
  generate_subFile(sim.dataType='linear-weak',
                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms
  generate_subFile(sim.dataType='nl-quadratic',
                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms with intercept
  generate_subFile(sim.dataType='nl-shiftquad',
                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Gene-gene interaction terms
  generate_subFile(sim.dataType='nl-interaction',
                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Sine interaction
  generate_subFile(sim.dataType='nl-sine',
                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
}


## ############################# BE10Asso00 ###############################

date = "071625"

for (norm_type in 0:6) {
  # Linear risk with moderate effects
  generate_subFile(sim.dataType='linear-moderate',
                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Linear risk with weak effects
  generate_subFile(sim.dataType='linear-weak',
                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms
  generate_subFile(sim.dataType='nl-quadratic',
                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms with intercept
  generate_subFile(sim.dataType='nl-shiftquad',
                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Gene-gene interaction terms
  generate_subFile(sim.dataType='nl-interaction',
                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Sine interaction
  generate_subFile(sim.dataType='nl-sine',
                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
}



## ############################# BE11Asso00 ###############################

date = "080425"

for (norm_type in 0:6) {
  # Linear risk with moderate effects
  generate_subFile(sim.dataType='linear-moderate',
                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Linear risk with weak effects
  generate_subFile(sim.dataType='linear-weak',
                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms
  generate_subFile(sim.dataType='nl-quadratic',
                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms with intercept
  generate_subFile(sim.dataType='nl-shiftquad',
                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Gene-gene interaction terms
  generate_subFile(sim.dataType='nl-interaction',
                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Sine interaction
  generate_subFile(sim.dataType='nl-sine',
                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
}


## ############################# BE10Asso10 ###############################

date = "072125"

for (norm_type in 0:6) {

  # Linear risk with moderate effects
  generate_subFile(sim.dataType='linear-moderate',
                   he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Linear risk with weak effects
  generate_subFile(sim.dataType='linear-weak',
                   he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms
  generate_subFile(sim.dataType='nl-quadratic',
                   he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms with intercept
  generate_subFile(sim.dataType='nl-shiftquad',
                   he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Gene-gene interaction terms
  generate_subFile(sim.dataType='nl-interaction',
                   he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Sine interaction
  generate_subFile(sim.dataType='nl-sine',
                   he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
}


## ############################# BE11Asso10 ###############################
date = "091625"

for (norm_type in 0:6) {
  
  # Linear risk with moderate effects
  generate_subFile(sim.dataType='linear-moderate',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
  
  # Linear risk with weak effects
  generate_subFile(sim.dataType='linear-weak',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
  
  # Squared terms
  generate_subFile(sim.dataType='nl-quadratic',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
  
  # Squared terms with intercept
  generate_subFile(sim.dataType='nl-shiftquad',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
  
  # Gene-gene interaction terms
  generate_subFile(sim.dataType='nl-interaction',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
  
  # Sine interaction
  generate_subFile(sim.dataType='nl-sine',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
                   date=date, GPUtype='l4')
}




## ############################# BE11Asso11 ###############################
date = "091725"

for (norm_type in 0:6) {

  # Linear risk with moderate effects
  generate_subFile(sim.dataType='linear-moderate',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Linear risk with weak effects
  generate_subFile(sim.dataType='linear-weak',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms
  generate_subFile(sim.dataType='nl-quadratic',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Squared terms with intercept
  generate_subFile(sim.dataType='nl-shiftquad',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Gene-gene interaction terms
  generate_subFile(sim.dataType='nl-interaction',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
                   date=date, GPUtype='l4')

  # Sine interaction
  generate_subFile(sim.dataType='nl-sine',
                   he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
                   date=date, GPUtype='l4')
}

