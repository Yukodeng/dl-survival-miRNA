# ########################## Run Simulation ##################################
# DATE     NOTE
# ############################################################################
options(echo=TRUE) # NOTE: supply job array number 1-7
idx <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# [NOTE] 3/5/26 for now leave out norm_type=2 (UQ)
## 4/6/26 for testing purpose also leave out norm type 3 (TMM) and 6 (Median)
norm_type_list <- c(0,1,3,4,5,6)
norm_type <- norm_type_list[idx]

dataType_vec <- c(
  "linear",
  "nl-quadratic",
  "nl-interaction",
  "nl-sine"
)
p_vec <- c(2, 5, 10, 30)
dataType_combs <- expand.grid(dataType_vec, p_vec, stringsAsFactors = FALSE)
dataType_list <- paste(dataType_combs$Var1, dataType_combs$Var2, sep = "-p")


# Source main pipeline script for raw data
source('~/dl-survival-miRNA/scripts/aug-batch-survival-pipeline-v3.R')

surv_folder = "sim_surv"
results_file = NULL #"model_results.csv"

# ## ###################### Test ###########################
# norm_type = 4
# purrr::map(dataType_list,
#   ~ sim.survdata(
#     sim.dataType = stringr::str_split(.x, "-p")[[1]][1],
#     p = stringr::str_split(.x, "-p")[[1]][2],
#     sim.data = sim.data,
#     train_sizes = c(100, 200, 500, 1000), #2000, 5000, 10000),
#     splits_per_size = c(3, 5, 5, 5), #10, 10, 10),
#     n_iter = 20,
#     he_train = 0,
#     he_test = 0,
#     beta_sort_train = 0,
#     beta_sort_test = 0,
#     norm_type = norm_type,
#     plot_km = T,
#     save_surv_data = F,
#     save_gene_data = F,
#     run_analysis = F,
#     stratify = 1,
#     save_results = T,
#     surv_folder = surv_folder,
#     results_file = results_file
#   )
# )

## ############################# BE00Asso00 ###############################
purrr::map(dataType_list,
  ~ sim.survdata(
    sim.dataType = stringr::str_split(.x, "-p")[[1]][1],
    p = stringr::str_split(.x, "-p")[[1]][2],
    sim.data = sim.data,
    train_sizes = c(100, 200, 500, 1000, 2000, 5000), #10000),
    splits_per_size = c(3, 5, 5, 5, 10, 10), #10),
    n_iter = 20,
    he_train = 0,
    he_test = 0,
    beta_sort_train = 0,
    beta_sort_test = 0,
    norm_type = norm_type,
    plot_km = F,
    save_surv_data = T,
    save_gene_data = T,
    run_analysis = F,
    stratify = 1,
    surv_folder = surv_folder,
    results_file = results_file
  )
)

## ############################# BE10Asso00 ###############################
purrr::map(dataType_list,
  ~ sim.survdata(
    sim.dataType = stringr::str_split(.x, "-p")[[1]][1],
    p = stringr::str_split(.x, "-p")[[1]][2],
    sim.data = sim.data,
    train_sizes = c(100, 200, 500, 1000, 2000, 5000), #10000),
    splits_per_size = c(3, 5, 5, 5, 10, 10), #10),
    n_iter = 20,
    he_train = 1,
    he_test = 0,
    beta_sort_train = 0,
    beta_sort_test = 0,
    norm_type = norm_type,
    plot_km = F,
    save_surv_data = T,
    save_gene_data = T,
    run_analysis = F,
    stratify = 1,
    surv_folder = surv_folder,
    results_file = results_file
  )
)

## ############################# BE11Asso00 ###############################
purrr::map(dataType_list,
  ~ sim.survdata(
    sim.dataType = stringr::str_split(.x, "-p")[[1]][1],
    p = stringr::str_split(.x, "-p")[[1]][2],
    sim.data = sim.data,
    train_sizes = c(100, 200, 500, 1000, 2000, 5000), #10000),
    splits_per_size = c(3, 5, 5, 5, 10, 10), #10),
    n_iter = 20,
    he_train = 1,
    he_test = 1,
    beta_sort_train = 0,
    beta_sort_test = 0,
    norm_type = norm_type,
    plot_km = F,
    save_surv_data = T,
    save_gene_data = T,
    run_analysis = F,
    stratify = 1,
    surv_folder = surv_folder,
    results_file = results_file
  )
)


## ################################# BE11Asso10 #####################################



## ################################# BE11Asso11 #####################################

