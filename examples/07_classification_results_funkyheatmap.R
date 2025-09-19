############################################################
##  Funky‑heatmap summary of Classification Accuracy results
##  ────────────────────────────────────────────────────────
##  • Pastel Morandi colours with **wider contrast range**
##    – gradients are now 6‑step and reach a darker anchor.
##  • Column‑group headers use an **even darker tint** of the
##    same hue, so every header bar is clearly visible.
##  • No dark greys / blacks; hues remain soft & JAMA‑friendly.
##  • Requires funkyheatmap ≥ 0.5.0  and  colorspace.
############################################################

rm(list = ls())

## ─── Packages ────────────────────────────────────────────
library(funkyheatmap)   # install.packages("funkyheatmap")
library(tidyverse)
library(colorspace)     # install.packages("colorspace")
library(abind)

## ─── Load pre‑computed result objects ───────────────────
load("./results/tables/classification_scenario1_analysis.RData")
clean_result = results
load("./results/tables/classification_scenario2_analysis.RData")
dirty_result = results

## ─── Helper: extract accuracy matrix, return long tibble ─
extract_accuracy_matrix <- function(result_list, scenario_name, param_name) {
  accuracies <- lapply(1:300, function(x) result_list[[param_name]][[x]]$test2$accuracy)
  array_data <- abind::abind(accuracies, along = 3)
  average_df <- apply(array_data, c(1, 2), mean)
  
  # Remove specific columns and update row/column names
  average_df2 <- t(average_df)[,-c(6,10)]
  rownames(average_df2) <- c("LASSO", "PAM", "KNN", "SVM", "SVM-FS", "RF", "RF-FS")
  colnames(average_df2) <- c("None", "TC", "UQ", "Med", "TMM", "DESeq", "PoissonSeq", 
                             "QN", "RUVg", "RUVr", "RUVs")

  # Convert to long format
  average_df2 |>
    as.data.frame() |>
    rownames_to_column("classification_method") |>
    pivot_longer(-classification_method,
                 names_to  = "normalization",
                 values_to = "accuracy") |>
    mutate(scenario   = scenario_name,
           parameters = param_name)
}

## ─── Extract representative parameter sets ──────────────
clean_c0.2     <- extract_accuracy_matrix(clean_result,
                                          "Biological",   "c0.2")
dirty_c0.2d1.5 <- extract_accuracy_matrix(dirty_result,
                                          "Real_world",   "c0.2_d1.5")

# Combine the available data
all_data <- bind_rows(clean_c0.2,
                      dirty_c0.2d1.5)

## ─── Prepare matrix for funkyheatmap ────────────────────
heatmap_data <- all_data |>
  unite("scenario_method", scenario, classification_method, sep = "_") |>
  select(normalization, scenario_method, accuracy) |>
  pivot_wider(names_from = scenario_method, values_from = accuracy) |>
  rename(id = normalization)

## ─── Pastel anchors and derived palettes ────────────────
# Base (light) pastel anchors for classification methods
anchor_hex <- c(
  "LASSO"  = "#F2B8B5",  # light rose
  "PAM"    = "#CFE4D0",  # mint pastel
  "KNN"    = "#D8D1E9",  # lavender pastel
  "SVM"    = "#F8DEC3",  # peach parchment
  "SVM-FS" = "#C7DFF3",  # powder blue
  "RF"     = "#E4D7C9",  # soft taupe
  "RF-FS"  = "#F5CFCB"   # salmon blush
)

## 1 | DARKER anchors for DATA cells (≈ 35 % darker)
cell_anchor <- colorspace::darken(anchor_hex, 0.35)

## 2 | Gradient palettes (white → darker anchor), 6 steps
cell_palettes <- purrr::map(cell_anchor,
                            ~ colorRampPalette(c("#FFFFFF", .x))(6)) |>
  set_names(tolower(names(cell_anchor)))

## 3 | EVEN DARKER anchors for HEADERS (≈ 50 % darker)
header_anchor <- colorspace::darken(anchor_hex, 0.50)
header_palettes <- purrr::map(header_anchor, \(x) rep(x, 2)) |>
  set_names(paste0(tolower(names(header_anchor)), "_hdr"))

## 4 | Combine all palettes
palettes <- c(cell_palettes, header_palettes)

## ─── Column & group metadata ────────────────────────────
classification_methods <- names(anchor_hex)
scenarios          <- c("Biological", "Real_world")
scenario_labels    <- c("Bio",        "Real")

### column_info: one row per heat‑map column
column_info <- tibble(
  id      = "id",
  group   = "normalization",
  name    = "",
  geom    = "text",
  palette = NA_character_,
  options = list(list())
)

for (method in classification_methods) {
  for (i in seq_along(scenarios)) {
    column_info <- add_row(
      column_info,
      id      = paste(scenarios[i], method, sep = "_"),
      group   = method,
      name    = scenario_labels[i],
      geom    = "funkyrect",
      palette = tolower(method),  # DATA cells use gradient palettes
      options = list(list())
    )
  }
}

### column_groups: header bars use *_hdr single‑colour palettes
column_groups <- tibble(
  group   = c("normalization", classification_methods),
  level1  = c("Normalization", classification_methods),
  palette = c(NA_character_,   paste0(tolower(classification_methods), "_hdr"))
)

## ─── Draw & save ────────────────────────────────────────
pdf("funky_heatmap_classification_summary.pdf", width = 16, height = 8)
funkyheatmap::funky_heatmap(
  data          = heatmap_data,
  column_info   = column_info,
  column_groups = column_groups,
  palettes      = palettes,
  scale_column  = FALSE,   # Accuracy already ranges 0‑1
  add_abc       = FALSE
)
dev.off()
