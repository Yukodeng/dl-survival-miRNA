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

library(funkyheatmap)
library(tidyverse)
library(colorspace)
library(abind)
library(ggplot2)
library(glue)

N_train = 5000

## ─── Load pre-computed result objects ───────────────────────

all_results <- read.csv(file.path("results", "all_results_w_stratified.csv"), check.name=F)
all_results <- all_results |> 
  mutate(
    `n train` = ifelse(
      (`model type` %in% c('rsf','stratified-rsf')) & (`n train`==8000), 10000, `n train`
    ),
    `model type` = gsub('svm',"ssvm",`model type`),
    `model type` = gsub('gb',"sgb",`model type`)
  ) 


# * Helper: return long tibble --------------------------- 

extract_cind_matrix <- function(all_results, batch_type) {

  results_sub  <- all_results |> 
    filter(`data type` != "nl-shiftquad") |> 
    mutate(
      Stratified = ifelse(grepl('stratified', `model type`), '_strat', ''),
      `model type` = gsub("stratified-", "", `model type`),
      Batch = gsub("_norm.*", "", `batchnorm type`),
      Normalization = gsub(".*_norm", "", `batchnorm type`),
      Normalization = case_when(
        Normalization == "None"     ~ "Raw",
        Normalization == "DEseq"    ~ "DESeq",
        Normalization == "Med"      ~ "Median",
        Normalization == "Quantile" ~ "QN",
        TRUE  ~ Normalization
      ),
      Association = case_when(
        `data type` == "linear-moderate" ~ "Moderate (linear)",
        `data type` == "linear-weak"     ~ "Weak (linear)",
        `data type` == "nl-quadratic"    ~ "Quadratic",
        # `data type` == "nl-shiftquad"    ~ "Centered Quadratic",
        `data type` == "nl-interaction"  ~ "Interactions", 
        `data type` == "nl-sine"         ~ "Sine"
      ),
      Model = case_when(
        `model type` == "deepsurv-torch" ~ "DeepSurv",
        `model type` == "oracle"         ~ "Oracle",
        `model type` == "oracle-linear"  ~ "Oracle (linear)",
        `model type` == "lasso"          ~ "LASSO",
        TRUE                             ~ toupper(`model type`)
      )
    ) |>
    filter( (`n train` == N_train) & (Batch == batch_type) ) |>
    group_by(Model, Stratified, Normalization, Association) |>
    summarize(`C-index` = mean(`test C`, na.rm = T)) |>
    ungroup() |>
    mutate(N = N_train,  Batch = batch_type)
      
  results_sub
}

batchTypes =  c(
  "BE00Asso00",
  "BE10Asso00",
  "BE11Asso00",
  "BE10Asso10",
  "BE11Asso10",
  "BE11Asso11"
)
result_list <- list()
for (batch_type in batchTypes) {
  result_list[[length(result_list)+1]] <- extract_cind_matrix(all_results, batch_type)
}

all_data <- do.call(bind_rows, result_list)

## fill in NA Oracle values with Oracle (linear) C-index for linear scenarios
for (norm in all_data$Normalization) {
  all_data[
    all_data$Model=='Oracle' & 
      grepl('linear', all_data$Association) &
      all_data$Normalization==norm, 
    "C-index"
  ] <- all_data |> 
        filter(Model=='Oracle (linear)' & grepl('linear',Association) & Normalization==norm) |> 
        pull(`C-index`)
}



## ─── Prepare matrix for funkyheatmap ────────────────────

# * Ordering -----------------------

association_types <- c("Moderate (linear)", "Weak (linear)","Quadratic", #"Centered Quadratic", 
                      "Interactions", "Sine")
association_labels <- association_types #LETTERS[1:6]
survival_methods <- c("Oracle", "Oracle (linear)", "LASSO", "DeepSurv", "RSF", "SSVM", "SGB")
normalizations <- c("Raw", "DESeq", "TMM", "TC", "Median", "UQ",  "QN")
norm_strat = c(normalizations, paste0(normalizations, "_strat"))
# norm_strat = unlist(Map(list, normalizations, paste0(normalizations, "_strat"))) |> as.vector()

# # * Transpose structure -----------

# all_wide <- all_data |>
#   unite("Normalization", Normalization, Stratified, sep = "") |>
#   mutate(
#     Batch         = factor(Batch, levels = batchTypes),
#     Association   = factor(Association, levels = association_types),
#     Normalization = factor(Normalization, levels = norm_strat),
#     Model         = factor(Model, levels = survival_methods)
#   ) |> 
#   unite("norm_assoc", Normalization, Association, sep = "_") |>
#   select(norm_assoc, Model, Batch, `C-index`) |>
#   pivot_wider(names_from = norm_assoc, values_from = `C-index`) |>
#   # unique row id per scenario × model
#   mutate(id = paste(Batch, Model, sep = " | "),  .after=Model) |> 
#   arrange(Batch, Model) #|> mutate(id = Model, .after=1) 



# # ─── Set color palettes ────────────────────────────────────

# anchor_hex <- c(
#   "Moderate (linear)"  = "#e9d4c4ff",
#   "Weak (linear)"      = "#eee2ddff",
#   "Quadratic"          = "#d8dbf4ff",
#   # "Centered Quadratic" = "#e1d1efff",
#   "Interactions"       = "#ffddc7ff",
#   "Sine"               = "#e5faddff"
#   # "Raw"    = "#a5a5a5ff",
#   # "TC"     = "#c6cefcff",
#   # "Median" = "#C7DFF3",
#   # "UQ"     = "#F2B8B5",
#   # "TMM"    = "#D8D1E9",
#   # "QN"     = "#CFE4D0",
#   # "DESeq"  = "#F8DEC3"
# )


# ## 1 | DARKER anchors for DATA cells (≈ 35 % darker)
# cell_anchor <- colorspace::darken(anchor_hex, 0.4)

# ## 2 | Gradient palettes (white → darker anchor), 6 steps
# cell_palettes <- purrr::map(cell_anchor,
#                             ~ colorRampPalette(c("#FFFFFF", .x))(6)) |>
#   set_names(tolower(names(cell_anchor)))

# ## 3 | EVEN DARKER anchors for HEADERS (≈ 40 % darker)
# header_anchor <- colorspace::darken(anchor_hex, 0.4)
# header_palettes <- purrr::map(header_anchor, \(x) rep(x, 2)) |>
#   set_names(paste0(tolower(names(header_anchor)), "_hdr"))

# ## 4 | Combine all palettes
# palettes <- c(cell_palettes, header_palettes)

# ## 5 | Add an overall palette for legend display
# palettes$funky_palette_grey <- RColorBrewer::brewer.pal(9, "Greys")[-1] |> rev()



# ## ─── Column Metadata: association × normalization type ────────────────────

# ### column_info: one row per heatmap column
# column_info <- tibble(
#   id      = "Model",#"id",
#   group   = "Association",
#   name    = "",
#   geom    = "text",
#   palette = NA_character_,
#   # ensure normalization IDs (left-most column) also use plain, 0° text
#   options = list(list(fontface = "plain", angle = 0))
# )

# for (i in seq_along(association_types)) {
#   for (norm in norm_strat) {#normalizations) {
#     column_info <- add_row(
#       column_info,
#       id      = paste(norm, association_types[i], sep = "_"),
#       group   = association_labels[i],
#       name    = norm,
#       geom    = "funkyrect",
#       palette = tolower(association_labels[i]),
#       # make the A/B/C/D plain, not italic, and non-angled
#       options = list(list(fontface = "plain", angle = 90))
#     )
#   }
# }

# ### column_groups: grouping by association type
# column_groups <- tibble(
#   Category   = c("Association", association_labels),
#   group     = c("Association", association_types),
#   palette   = c(NA_character_, paste0(tolower(association_labels), "_hdr"))
# )



# ## ─── Row Metadata: Model and batch type ─────────────────────────────--

# row_info <- all_wide |>
#   transmute(
#     id    = id, #as.character(Model),
#     group = as.character(Batch),  # scenario label for grouping
#     name  = Model,     # show only the model in the row label
#     geom  = "text"
#   )

# row_groups <- tibble(
#   group   = levels(all_wide$Batch),
#   level1  = levels(all_wide$Batch)#,
#   # palette = rep("scenario_hdr", length(levels(all_wide$batchnorm_type)))
# )


# ## ─── Legend  ────────────────────────────────────────────────

# ### Hide auto-generated legends for column groups and display the overall
# legends <- list(
#   list(
#     palette = "moderate (linear)",
#     enabled = F
#   ),
#   list(
#     palette = "weak (linear)",
#     enabled = F
#   ),
#   list(
#     palette = "quadratic",
#     enabled = F
#   ),
#   # list(
#   #   palette = "centered quadratic",
#   #   enabled = F
#   # ),
#   list(
#     palette = "interactions",
#     enabled = F
#   ),
#   list(
#     palette = "sine",
#     enabled = F
#   ),
#   list(
#     palette = "funky_palette_grey",
#     geom = "funkyrect",
#     title = "C-index",
#     enabled = TRUE,
#     labels = c("0.5", "", "0.6", "", "0.7", "", "0.8", "", "0.9", "", "1")
#   )
# )

# ### select "Model" column specified as id in column_info to display as row names
# heatmap_data <- all_wide |> 
#   mutate(Model=as.character(Model), .after = 1) |> 
#   select(-c(Batch))

# ### scale values to fit with range (0.5, 1)
# for (col in colnames(heatmap_data)[!colnames(heatmap_data) %in% c('id',"Model","Batch")]) {
#   heatmap_data[col] <- pmax(0, (heatmap_data[[col]] - 0.5)/0.5)
# }



# ## ─── Draw plot  ────────────────────────────────────────────────

# # pdf("funky_heatmap_summary_reoriented.pdf", width = 12, height = 15)
# svg(
#   here::here("results","plots", glue::glue("funky_heatmap_survival_results_wStrat_{N_train}.svg")), 
#   width = 15, height = 12
# )
# p <- funky_heatmap(
#   data          = heatmap_data,
#   column_info   = column_info,
#   column_groups = column_groups,
#   row_info      = row_info,
#   row_groups    = row_groups,
#   palettes      = palettes,
#   legends       = legends,
#   position_args = position_arguments(col_annot_angle = 60),
#   scale_column  = FALSE,
#   add_abc       = FALSE
# )
# p <- p + theme(
#   text = element_text(size=15, face = "plain"),
#   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
#   axis.text.y = element_text(face = "bold")
# )
# print(p)
# dev.off()




## ---- Test Stratified -----------------------------------------------

# prep_wide <- function(df, tag){
#    df |>  
#     mutate(
#       Association   = factor(Association, levels = association_types),
#       Normalization = factor(Normalization, levels = normalizations),
#       Model         = factor(Model, levels = survival_methods)
#     ) |>
#     unite("norm_assoc", Normalization, Association, sep="_") |>
#     mutate(norm_assoc = paste0(norm_assoc, "_", tag)) |>
#     select(Batch, Model, norm_assoc, `C-index`) |> View()
#     pivot_wider(names_from = norm_assoc, values_from = `C-index`)
# }

# wide_ns <- prep_wide(all_data |> filter(Stratified==""),       "NS")
# wide_st <- prep_wide(all_data |> filter(Stratified=="_strat"), "ST")

# Join NS/ST on Batch, Model and build unique row id + display label
all_wide <- all_data |> 
    mutate(
      Association   = factor(Association, levels = association_types),
      Normalization = factor(paste0(Normalization, Stratified), levels = norm_strat),
      Model         = factor(Model, levels = survival_methods)
    ) |>
    unite("norm_assoc", Normalization, Association, sep="_") |>
    select(Batch, Model, norm_assoc, `C-index`) |> 
    pivot_wider(names_from = norm_assoc, values_from = `C-index`) |> 
    # full_join(wide_ns, wide_st, by = c("Batch","Model")) |> 
    arrange(Batch, Model) |> 
    mutate(
      id = paste(Batch, Model, sep = " | "), .after=Model
    ) 



# ================== 2) Column metadata (NS circle, ST square) ==================
# Left label column
column_info <- tibble(
  id      = "Model",
  group   = "Association",
  name    = "",
  geom    = "text",
  palette = NA_character_
)

# Build two columns under each Association × Normalization: ..._NS and ..._ST
for (a in association_types) {
  for (n in norm_strat) {
    base <- paste(n, a, sep = "_")
    # NS
    column_info <- add_row(    # <-- non-stratified
      column_info,
      id       = base,
      group    = a,
      name     = gsub("_strat", "", n),
      geom     = "funkyrect",     
      palette  = tolower(a)
    )
    # # ST
    # column_info <- add_row(    # <-- stratified
    #   column_info,
    #   id       = paste0(base, "_ST"),
    #   group    = a,
    #   name     = n,  #paste0(n, "-ST"),
    #   geom     = "funkyrect",       
    #   palette  = tolower(a)
    # )
  }
}

column_groups <- tibble(
  group   = c("Association", association_types),
  level1  = c("Association", association_labels),
  palette = c(NA_character_, paste0(tolower(association_types), "_hdr"))
)



# ==========  Legend ==================================
# Association headers (unchanged)
anchor_hex <- c(
  "Moderate (linear)"  = "#e1d0c7",
  "Weak (linear)"      = "#dfdfdf",
  "Quadratic"          = "#ced6f9",
  "Interactions"       = "#f7cfae",
  "Sine"               = "#CFE4D0"
)

## 1 | DARKER anchors for DATA cells (≈ 35 % darker)
cell_anchor <- colorspace::darken(anchor_hex, 0.4)

## 2 | Gradient palettes (white → darker anchor), 6 steps
cell_palettes <- purrr::map(cell_anchor,
                            ~ colorRampPalette(c("#FFFFFF", .x))(6)) |>
  set_names(tolower(names(cell_anchor)))

## 3 | EVEN DARKER anchors for HEADERS (≈ 40 % darker)
header_anchor <- colorspace::darken(anchor_hex, 0.4)
header_palettes <- purrr::map(header_anchor, \(x) rep(x, 2)) |>
  set_names(paste0(tolower(names(header_anchor)), "_hdr"))

## 4 | Combine all palettes
palettes <- c(cell_palettes, header_palettes)

## 5 | Add an overall palette for legend display
palettes$funky_palette_grey <- RColorBrewer::brewer.pal(9, "Greys")[-1] |> rev()

## 5 | Add legend for stratification
palettes <- c(
  palettes,
  # a tiny categorical palette just to create a "Stratification" legend (shape cue)
  list(strata_pal = c(NS = "#B0B0B0", ST = "#B0B0B0"))
)



# ================== 3) Row grouping (batches) ==================
row_info <- all_wide |>
  transmute(
    id = id, 
    group = as.character(Batch),
    name  = Model,     # show only the model in the row label
    geom  = "text"
  )

row_groups <- tibble(
  group   = levels(factor(all_wide$Batch)),
  level1  = levels(factor(all_wide$Batch))
)


# ================== 4) Legends ================================
# C-index legend (continuous grey) + a shape legend for stratification
legends <-list(
  list(
    palette = "moderate (linear)",
    enabled = F
  ),
  list(
    palette = "weak (linear)",
    enabled = F
  ),
  list(
    palette = "quadratic",
    enabled = F
  ),
  list(
    palette = "interactions",
    enabled = F
  ),
  list(
    palette = "sine",
    enabled = F
  ),
  list(
    palette = "funky_palette_grey",
    geom = "funkyrect",
    title = "C-index",
    enabled = TRUE,
    labels = c("0.5", "", "0.6", "", "0.7", "", "0.8", "", "0.9", "", "1")
  )#,

  # # two small “dummy” entries to show shapes for NS vs ST
  # list(
  #   palette = "strata_pal",
  #   title = "Stratification",
  #   geom = "funkyrect",
  #   labels = "Non-stratified"
  # ),
  # list(
  #   palette = "strata_pal", 
  #   title = "Stratification",
  #   geom = "funkyrect",
  #   labels = "Stratified"
  # )
)# |> unlist(recursive = FALSE)




# Make size columns (0.5 -> 0; 1.0 -> 1)
heatmap_data <- all_wide |> 
  mutate(Model=as.character(Model)) |> 
  select(-c(Batch))

### scale values to fit with range (0.5, 1)
for (col in colnames(heatmap_data)[!colnames(heatmap_data) %in% c('id',"Model","Batch")]) {
  heatmap_data[col] <- pmax(0, (heatmap_data[[col]] - 0.5)/0.5)
}


# ================== 5) Draw ==================


p <- funkyheatmap::funky_heatmap(
  data          = heatmap_data,
  column_info   = column_info,
  column_groups = column_groups,
  row_info      = row_info,
  row_groups    = row_groups,
  palettes      = palettes,
  legends       = legends,
  scale_column  = FALSE,
  add_abc       = FALSE,
  position_args = position_arguments(
    col_annot_angle = 60,
    col_annot_offset = 4,
    expand_xmax = 4
  ) 
)

p <- p + theme(
  text         = element_text(size = 7),
  axis.text.x  = element_text(size = 6, hjust = 1),
  axis.text.y  = element_text(size = 7, face = "bold")
)

svg(
  here::here("examples", glue("example_funkyheatmap_wStrat_{N_train}.svg")), 
  width = 15, height = 12
)
print(p)
dev.off()
