############################################################
##  Funky-heatmap summary of C-index results 
############################################################
# DATE       NOTES
# 08OCT2025  Created script from script shared by Jian
# 09OCT2025  A few changes: 
#            1) remove nl-shiftquad condition
#            2) change ordering of normalization methods (raw,DESeq,TMM,TC,Med,UQ/QN)
#            3) relabel models as SSVM / SGB
#            4) figure out display of stratified models (Overlay?)
############################################################

rm(list = ls())

## ─── Packages ────────────────────────────────────────────

# library(funkyheatmap)
library(tidyverse)
library(colorspace)
library(abind)
library(ggplot2)
library(ggh4x)
library(glue)


## ---- organize model results ----

setwd("~/dl-survival-miRNA")
root_dir <- file.path('results', 'models')

batchTypes = c("BE00Asso00", "BE10Asso00", "BE11Asso00")
normTypes = c("None")#, "TC", "Quantile", 'DEseq')#, "TMM", "Quantile", "DEseq", "Med") 
batchNorm_combs <- expand.grid(batchTypes, normTypes, stringsAsFactors = FALSE)
batchType_folders <- paste(batchNorm_combs$Var1, batchNorm_combs$Var2, sep = "_norm")

dataType_vec <- c(
  "linear",
  "nl-quadratic",
  "nl-interaction",
  "nl-sine"
)
p_vec <- c(5, 10, 30)#, 2, 60)
dataType_combs <- expand.grid(dataType_vec, p_vec, stringsAsFactors = FALSE)
dataType_folders <- paste(dataType_combs$Var1, dataType_combs$Var2, sep = "-p")
dataType_folders <- c(dataType_folders, 
  c("linear-p5-2x",
    "nl-quadratic-p5-2x",
    "nl-interaction-p5-2x",
    "nl-sine-p5-2x"),
  c("linear-p30_2",
    "nl-quadratic-p30_2",
    "nl-interaction-p30_2",
    "nl-sine-p30_2")
)

  
# ---- find all CSVs under */cox_all/* within the selected batch folders ----

split_batch_norm_cols <- function(
  df, 
  key_col = "batch_type", 
  from_col = "norm_type", 
  to_cols = c("batch_type", "norm_type"), sep = "_") {
  
  if (!key_col %in% colnames(df)) {
    df |> 
      separate(get(from_col), into = to_cols, sep = sep, remove = T)
  } else {
    df
  }
}


load_results_from_path <- function(batchType, dataType, modelType, pattern="all", cols_rm, sep = "\t") {

  path = file.path(root_dir, batchType, dataType, modelType)
  f = list.files(path = path, pattern = pattern, full.names = T)
  f = f[!grepl("archive_*", f)]

  if (length(f)==0) {
    message(glue::glue("File not found in folder {path}"))
    NULL
  } else {
    f = f[length(f)]
    message(f)
    df <- read.csv(f, sep = sep)
    if (ncol(df)==1) {
      df <- read.csv(f, sep = setdiff(c("\t", ","), sep))
    }
    df <- df %>% select(all_of(colnames(.)[!colnames(.) %in% cols_rm]))

    split_batch_norm_cols(df)
  }

}

# ---- read all results ----

all_results_list <- list()

for (batchType in batchType_folders) {
  for (dataType in dataType_folders) {

    cox_df <- load_results_from_path(batchType, dataType, 'cox', cols_rm = c("event_rate", "lambda"))
    cox_2x_df <- load_results_from_path(batchType, dataType, 'cox', pattern = "p5_2x", cols_rm = c("event_rate", "lambda"))
    cox_30_2_df <- load_results_from_path(batchType, dataType, 'cox', pattern = "p30_v2", cols_rm = c("event_rate", "lambda"))

    # ds_df <- load_results_from_path(batchType, dataType, 'deepsurv-torch', cols_rm="train_time", sep=",")
    # ds_app_df <- load_results_from_path(batchType, dataType, 'deepsurv-torch', pattern="042126", cols_rm="train_time", sep=",")
    # # ds_strat_df <- load_results_from_path(batchType, dataType, 'stratified-deepsurv-torch', "train_time", ",")
    # sgb_df <- load_results_from_path(batchType, dataType, 'sgb', cols_rm="train_time", sep=",")
    # sgb_strat_df <- load_results_from_path(batchType, dataType, 'stratified-sgb', cols_rm="train_time", sep=",")

    if (all(sapply(
      list(cox_df, cox_2x_df, cox_30_2_df), #ds_df, ds_app_df, sgb_df, sgb_strat_df), 
      is.null))) next

    result_df <- rbind(
      cox_df, cox_2x_df, cox_30_2_df
      # ds_df, ds_app_df,
      # sgb_df, sgb_strat_df
    )
    
    if (length(result_df) != 0) {
      all_results_list[[length(all_results_list)+1]] <- result_df |> data.frame()
    }
  }
}

# ---- stack ----
all_results <- do.call(rbind, all_results_list)


# # ----- [TEMPORARY] Test jobs results retrieval -----
# results_df = load_results_from_path('BE11Asso00_normNone', "linear-weak", 'cox', "5x", c("event_rate", "lambda"))
# results_df |> 
#   filter(dataset == 'test', stratified == TRUE) |> 
#   group_by(
#     batch_type, norm_type, data_type, stratified, model, n_train
#   ) |> 
#   summarize(`C-index` = mean(cind, na.rm = T),
#             min = min(cind, na.rm = T),
#             max = max(cind, na.rm = T)) |> 
#   ungroup() |> View()



## ─── Prepare Dataset  ───────────────────────────────────────


library(stringr)

# n_train_sel = 5000
# batch_sel = "BE10Asso00"
# modelType_sel = "lasso" # "oracle" "oracle (linear)" "lasso" "sgb"
dataType_sel = "sine" # "linear" "quadratic" "interaction" "sine"


data = all_results |> 
  filter(dataset == 'test',
        grepl(dataType_sel, data_type)
  ) |>
  mutate(model = gsub("stratified-", "", model)) |>
  mutate(
    Stratified = ifelse(stratified,'Stratified', 'Non-stratified'),
    Model = factor(model,
      levels = c("oracle-nl", "oracle-linear", "lasso"), #"sgb", "deepsurv-torch"),
      labels = c("oracle", "oracle (linear)", "lasso"), #"sgb", "deepsurv")
    ),
    Batch = factor(batch_type,
      levels = c('BE00Asso00', 'BE10Asso00', 'BE11Asso00'),
      labels = c('BE00', 'BE10', 'BE11')),
    Normalization = gsub("norm", "", norm_type),
    Normalization = case_when(
      Normalization == "None"     ~ "Raw",
      Normalization == "DEseq"    ~ "DESeq",
      Normalization == "Med"      ~ "Median",
      Normalization == "Quantile" ~ "QN",
      TRUE  ~ Normalization
    ),
    Normalization = factor(
      Normalization#, 
      # levels = c("Raw", "TC", "QN", "DESeq")#, "TMM", "TC", "Median", "QN")
    ),
    Association = factor(case_when(
        # data_type == "linear-p2"    ~ "Linear (p=2)",
        data_type == "linear-p5-2x" ~ "Linear (p=5; 2x)",
        data_type == "linear-p5"    ~ "Linear (p=5)",
        data_type == "linear-p10"   ~ "Linear (p=10)",
        data_type == "linear-p30"   ~ "Linear (p=30)",
        data_type == "linear-p30_2" ~ "Linear (p=30; v2)",
        # data_type == "linear-p60" ~ "Linear (p=60)",

        # data_type == "nl-quadratic-p2"     ~ "Quadratic (p=2)",
        data_type == "nl-quadratic-p5-2x"  ~ "Quadratic (p=5; 2x)",
        data_type == "nl-quadratic-p5"     ~ "Quadratic (p=5)",
        data_type == "nl-quadratic-p10"    ~ "Quadratic (p=10)",
        data_type == "nl-quadratic-p30"    ~ "Quadratic (p=30)",
        data_type == "nl-quadratic-p30_2"  ~ "Quadratic (p=30; v2)",
        # data_type == "nl-quadratic-p60" ~ "Quadratic (p=60)"

        # data_type == "nl-interaction-p2"  ~ "Interaction (p=2)",
        data_type == "nl-interaction-p5-2x" ~ "Interaction (p=5; 2x)",
        data_type == "nl-interaction-p5"  ~ "Interaction (p=5)",
        data_type == "nl-interaction-p10" ~ "Interaction (p=10)",
        data_type == "nl-interaction-p30" ~ "Interaction (p=30)",
        data_type == "nl-interaction-p30_2" ~ "Interaction (p=30; v2)",

        # data_type == "nl-sine-p2"    ~ "Sine (p=2)",
        data_type == "nl-sine-p5-2x" ~ "Sine (p=5; 2x)",
        data_type == "nl-sine-p5"    ~ "Sine (p=5)",
        data_type == "nl-sine-p10"   ~ "Sine (p=10)",
        data_type == "nl-sine-p30"   ~ "Sine (p=30)",
        data_type == "nl-sine-p30_2" ~ "Sine (p=30; v2)"
      ),
      levels = c(
        "Linear (p=5; 2x)", "Linear (p=5)", "Linear (p=10)", "Linear (p=30)", "Linear (p=30; v2)", 
        "Quadratic (p=5; 2x)", "Quadratic (p=5)", "Quadratic (p=10)","Quadratic (p=30)","Quadratic (p=30; v2)",
        "Interaction (p=5; 2x)","Interaction (p=5)", "Interaction (p=10)","Interaction (p=30)","Interaction (p=30; v2)",
        "Sine (p=5; 2x)", "Sine (p=5)", "Sine (p=10)", "Sine (p=30)", "Sine (p=30; v2)"
      )
    )
  ) |>
  group_by(
    Batch, Normalization,  Association, Stratified, Model, n_train#`n train`
  ) |> 
  summarize(`C-index` = mean(cind, na.rm = T),
            min = min(cind, na.rm = T),
            max = max(cind, na.rm = T)) |> 
  ungroup() |>
  select(
    N = n_train,
    Model,
    Batch,
    Normalization, 
    Stratified,
    Association, 
    `C-index`, #sd,
    min, max
  )



## ─── Draw Bar Plot ──────────────────────────────────────────

pd <- position_dodge(width = 1)
p <- data |>
  filter(
    !N %in% c(1000,2000,5000)
    # as.character(Model) == modelType_sel
  ) |> 
  ggplot(
    aes(x = Model, y = `C-index`, color = Model, group = interaction(Model, Stratified))) +
  geom_errorbar(
    aes(ymin = min, ymax = max, linetype=Stratified), 
    color = "grey",
    width = 0.3, linewidth = .75, position = pd
  ) +
  geom_point(aes(shape = Stratified), size = 2.75, position = pd) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color='grey20', alpha=0.3)+
  
  facet_nested(Batch ~ Association + N) +
  # facet_grid(rows = vars(Association), cols = vars(Batch), space = "free_x") +

  scale_color_brewer(palette = "Set2") + 
  scale_shape_manual(values = c(`Non-stratified` = 16, `Stratified` = 25)) +
  scale_linetype_manual(values = c(`Non-stratified` = "solid", `Stratified` = "solid")) +
  labs(
    x = NULL,
    y="Test C-index", 
    title = glue("Simulation results for {dataType_sel} associations (No Normalization)") #glue("{modelType} Results")
  ) + 
  # scale_y_continuous(breaks = seq(0.5, .9, by=0.1)) + 
  ylim(0.45, 0.95) +
  theme_bw(base_size = 18, base_line_size = .3, base_rect_size = .2) +
  theme(
    panel.spacing = unit(0, "lines"),
    strip.text.y = element_text(size = 13, angle = 0),
    axis.text.x = element_blank(), #element_text(size=10, angle = 55, hjust=1),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.margin = margin(t=-5)
  ) + 
  guides(
    color    = guide_legend(title = "Model", nrow = 1, order = 1),
    shape    = guide_legend(title = "Stratified", nrow = 1, order = 2),
    linetype = guide_legend(title = "Stratified", nrow = 1, order = 2)
  )
p

## ─── Save Plot ───────────────────────────────────────────

ggsave(
  file.path("results", "plots", "2026", 
            # glue("facet_scatterplot_{dataType_sel}_normNone_results_v3.jpg")
            glue("test_30_v2_scatterplot_{dataType_sel}_normNone_results_v3.jpg")
          ), 
  plot=p, height=8, width=16, units='in', dpi=600
)
