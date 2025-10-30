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


modelType = "sgb" # "oracle", "oracle-linear", "lasso", "deepsurv-torch", "rsf", "ssvm", "sgb"



## ─── Load pre-computed result objects ───────────────────────

all_results <- read.csv(file.path("results", "all_results_w_stratified.csv"), check.name=F)
all_results <- all_results |> 
  mutate(
    # `n train` = ifelse(
    #   (`model type` %in% c('rsf','stratified-rsf')) & (`n train`==8000), 10000, `n train`
    # ),
    `model type` = gsub('svm',"ssvm",`model type`),
    `model type` = gsub('gb',"sgb",`model type`)
  ) 

if (modelType=='rsf') {
  all_results <- all_results |>  filter(`n train` != 10000)
}


## fill in NA Oracle values with Oracle (linear) C-index for linear scenarios
for (norm in unique(all_results$`batchnorm type`)) {
  all_results[
      all_results$`model type` %in% c('oracle', 'stratified-oracle') & 
      grepl('linear', all_results$`data type`) &
      all_results$`batchnorm type`==norm, 
    c("train C", "test C")
  ] <- 
    all_results |> 
      filter(
        `model type` %in% c('oracle-linear', 'stratified-oracle-linear') & 
        grepl('linear', `data type`) & 
        `batchnorm type`==norm) |> 
      select(`train C`, `test C`)
}



## ─── Prepare Dataset  ───────────────────────────────────────

data = all_results |> 
  filter(`data type` != "nl-shiftquad") |> 
  mutate(
    Stratified = factor(
      ifelse(grepl('stratified', `model type`), 'ST', 'NS')
    ),
    Model = gsub("stratified-", "", `model type`),
    Batch = gsub("_norm.*", "", `batchnorm type`),
    Normalization = gsub(".*_norm", "", `batchnorm type`),
    Normalization = case_when(
      Normalization == "None"     ~ "Raw",
      Normalization == "DEseq"    ~ "DESeq",
      Normalization == "Med"      ~ "Median",
      Normalization == "Quantile" ~ "QN",
      TRUE  ~ Normalization
    ),
    Normalization = factor(
      Normalization, 
      levels = c("Raw", "DESeq", "TMM", "TC", "Median", "UQ",  "QN")
    ),
    Association = factor(case_when(
        `data type` == "linear-moderate" ~ "Moderate (linear)",
        `data type` == "linear-weak"     ~ "Weak (linear)",
        `data type` == "nl-quadratic"    ~ "Quadratic",
        `data type` == "nl-interaction"  ~ "Interactions", 
        `data type` == "nl-sine"         ~ "Sine"
      ),
      levels = c("Moderate (linear)","Weak (linear)","Quadratic","Interactions","Sine")
    )
  ) |>
  filter(#Batch == batchType & 
    Model == modelType) |>
  group_by(Batch, Normalization, `n train`, Association, Stratified) |> 
  summarize(`C-index` = mean(`test C`, na.rm = T),
            min = min(`test C`, na.rm = T),
            max = max(`test C`, na.rm = T)) |> 
  ungroup() |> 
  select(
    N = `n train`,
    Batch,
    Normalization, 
    Stratified,
    Association, 
    `C-index`, #sd,
    min, max
  )




## ─── Draw Bar Plot ──────────────────────────────────────────

pd <- position_dodge2(width = .4 ,preserve = "single")
p <- data |>
  ggplot(
    aes(x = Normalization, y = `C-index`, color = Normalization, group = interaction(Normalization, Stratified))) +
  geom_errorbar(
    aes(ymin = min, ymax = max, linetype=Stratified), 
    color = "grey",
    width = 0.3, linewidth = .75, position = pd
  ) +
  geom_point(aes(shape = Stratified), size = 1.75, position = pd) +
  
  facet_nested(Batch ~ Association + N) +
  # facet_grid(rows = vars(Association), cols = vars(N), space = "free_x") +

  scale_color_brewer(palette = "Set2") + 
  scale_shape_manual(values = c(NS = 16, ST = 25)) +
  scale_linetype_manual(values = c(NS = "solid", ST = "dashed")) +
  labs(
    x = NULL,
    y="C-Index", 
    title = glue("{modelType} Results")
  ) + 
  # scale_y_continuous(breaks = seq(0.4, 1, by=0.2)) + #ylim(0.45, 0.95) +
  theme_bw(base_size = 14, base_line_size = .3, base_rect_size = .2) +
  theme(
    panel.spacing = unit(0, "lines"),
    strip.text.y = element_text(size = 10, angle = 0),
    axis.text.x = element_blank(), #element_text(size=10, angle = 55, hjust=1),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
  ) + 
  guides(
    color    = guide_legend(title = "Normalization", nrow = 1, order = 1),
    shape    = guide_legend(title = "Stratified", nrow = 1, order = 2),
    linetype = guide_legend(title = "Stratified", nrow = 1, order = 2)
  )


## ─── Save Plot ───────────────────────────────────────────

ggsave(
  file.path("results", "plots", glue("facet_boxplot_{modelType}_results.jpg")), 
  plot=p, height=9, width=20, units = 'in'
)
