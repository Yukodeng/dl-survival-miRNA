############################################################################# #
# Date      Note
# 24FEB2026 Run for new augmented data (300000) with 15 batches and partial constant BE
# 05MAR2026 Move survival simulation to separate script and generate all at once per asso type
#           Change to partition samples per train size and iter id to ensure non-overlap
#           Save results per experiment condition per iteration by appending to file 
############################################################################# #
library(ggsurvfit)
library(survival)
library(glmnet)
library(edgeR)
library(dplyr)
library(pec)
library(ggplot2)
library(cli)
library(nanoparquet)
library(preprocessCore)
options(preprocessCore.nthreads = 1)
options(stringsAsFactors = F)
Sys.setenv("RCPP_PARALLEL_BACKEND" = "single")

# Set working directory and load functions ----------------------------------

setwd("~/dl-survival-miRNA/")
source("scripts/utils/norm.R")
# [NOTE] First time: make sure to run script for simulating survial outcome
# source("scripts/simulate-survival-outcome.R")

set.seed(1234)

# Load Raw Augmented Gene Counts --------------------------------------------

g.names = read.csv("data/augmentation_miRNA_names_538.csv")$gene # gene names

sim.data = read_parquet(file.path(
  "data", "simGeneExp_full_preprocessed.parquet"
))

# ---- Helper functions ----------------------------------------------------

## --- Normalization map and implement tool ----
norm_map <- list(
  `0` = list(name = "None", fun = NULL),
  `1` = list(name = "TC", fun = norm.TC),
  `2` = list(name = "UQ", fun = norm.UQ),
  `3` = list(name = "TMM", fun = norm.TMM),
  `4` = list(name = "Quantile", fun = NULL),
  `5` = list(name = "DEseq", fun = norm.DESeq),
  `6` = list(name = "Med", fun = norm.med)
)

apply_sf_then_log <- function(x_count, sf) {
  x_count_norm <- t(t(x_count) / sf)
  log(x_count_norm + 0.5)
}

## ---- Train and test set indexing -----------
alloc_counts <- function(n, n_batch, iter_i) {
  # Split n_total across batches: base or base+1, rotating which batches get +1 by iter_i
  batches <- seq_len(n_batch)
  base <- n %/% n_batch
  r <- n %% n_batch
  start <- (iter_i %% n_batch) + 1
  order <- c(batches[start:n_batch], batches[1:(start - 1)])
  plus_one <- order[seq_len(r)]

  setNames(base + as.integer(batches %in% plus_one), batches)
}

# train_size=500
# test_size=1000
# iter_i=21
# n_batch = 15
# max_train_per_batch = ceiling(5000 / n_batch)
# N = 20000

extract_split_ids <- function(train_size, test_size, iter_i, max_train_per_batch, n_batch = 15, N = 20000) {
  key <- paste0("train", train_size, "_iter", iter_i)

  max_test_per_batch <- ceiling(test_size / n_batch)
  chunk_len <- max_train_per_batch + max_test_per_batch
  train_counts <- alloc_counts(train_size, n_batch, iter_i)
  test_counts <- alloc_counts(test_size, n_batch, iter_i)

  start <- (iter_i - 1) * chunk_len + 1
  end <- iter_i * chunk_len
  test_start <- max_train_per_batch + 1
  train_idx <- integer(0)
  test_idx <- integer(0)

  for (b in 1:n_batch) {
    chunk <- ((b - 1) * N + start):((b - 1) * N + end)
    tr_take <- unname(train_counts[b])
    te_take <- unname(test_counts[b])
    
    # train: from front of chunk
    train_idx_b <- chunk[seq_len(tr_take)]
    # test: from fixed idx to avoid any overlap with train across train sizes
    test_idx_b <- chunk[test_start:(test_start + te_take - 1)]

    train_idx <- c(train_idx, train_idx_b)
    test_idx <- c(test_idx, test_idx_b)
  }
  list(
    key = key,
    train_idx = train_idx,
    test_idx = test_idx,
    train_counts = train_counts,
    test_counts = test_counts
  )
}

## ---- Transform dat to nonlinear form -----------
x0_transform <- function(sim.dataType, x0, interact.gene = NULL) {
  if (grepl("linear", sim.dataType)) { #%in% c('linear-moderate', 'linear-weak')) {
    return(NULL)
  }
  ## 10 quadratic terms -------------------
  if (sim.dataType == 'nl-quadratic') {
    x0_t = apply(x0, 2, scale)^2 ## scale across 2=columns (genes)

    ## 10 quadratic with intercept ----------
  } else if (sim.dataType == 'nl-shiftquad') {
    x0_scaled = apply(x0, 2, scale)
    c = apply(x0_scaled, 2, median)
    x0_t = sweep(x0_scaled, 2, c, "-")^2 # element-wise x0_tansformation

    ## 10 interaction terms ----------------
  } else if (sim.dataType == 'nl-interaction') {
    x0_scaled = apply(x0, 2, scale)
    interact.gene.x = x0_scaled[, interact.gene]
    x0_t = sweep(x0_scaled, 1, interact.gene.x, FUN = "*")

    ## 10 sine terms ----------------------
  } else if (sim.dataType == 'nl-sine') {
    x0_t = sin(x0)
  } else {
    stop("Unknown sim.dataType: ", sim.dataType)
  }
  return(as.data.frame(x0_t))
}

## ----- plot KM curves -----
plot_KM <- function(df) {
  fig <- survfit(
    Surv(time, status) ~ type,
    data = df,
    conf.type = "log-log") |>
    ggsurvfit() +
    labs(
      title = 'Survival',
      x = 'Time',
      y = 'Survival Probability'
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    add_confidence_interval() +
    add_risktable() +
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    theme_bw()

  print(fig)
}

## ----- Fit CoxPH with ridge fallback and apply to test ----
# X=x0_train
# y=y_train
# beta_init = NULL
# theta = 1
# stratify = T
fit_coxph_safe <- function(formula, X, y, beta_init = NULL, theta = 1, stratify = 0) {
  df <- data.frame(cbind(X, y))
  X <- as.matrix(X)
  if (!is.null(beta_init)) {
    fit <- coxph(
      formula = formula,
      data = df,
      x = TRUE,
      y = TRUE,
      init = beta_init,
      control = coxph.control(iter.max = 0)
    )
  } else {
    fit <- tryCatch(
      coxph(formula = formula, data = df, x = TRUE, y = TRUE),
      warning = function(w) w,
      error = function(e) e
    )
    if (inherits(fit, "warning") || inherits(fit, "error")) {
      # rhs <- paste0( "ridge(", paste(colnames(X), collapse = ","), ", theta=", theta,")")
      rhs <- paste0( "ridge(X, theta=", theta,")")
      if (stratify) {
        rhs <- paste(rhs, "+ strata(batch_id)")
      }
      f2 <- as.formula(paste("Surv(time, status) ~", rhs))
      fit <- coxph(formula = f2, data = df, x = TRUE, y = TRUE)
    }
  }
  fit
}

## ----- Evaluation functions ------------------------
get_IBS_score <- function(
  df,
  fit,
  time_col = "time",
  status_col = "status",
  k_times = 20
) {
  stopifnot(is.data.frame(df))
  # stopifnot(all(
  #   all.vars(fit$formula) %in% names(df) | all.vars(fit$formula) %in% c("Surv")
  # ))
  # Time grids: k_times equally spaced points between the min and max EVENT times
  events <- df[[status_col]] == 1
  times <- seq(
    ceiling(min(df[events, time_col])),
    floor(max(df[events, time_col])),
    length.out = k_times
  )
  # Compute Brier curves and IBS (IPCW with KM for censoring)
  ibs <- tryCatch({
    pec <- pec(
        object = list('coxph' = fit),
        formula = f.pec, #fit$formula,
        data = df,
        times = times,
        exact = F,
        reference = F,
        cens.model = "cox",
        keep.matrix = F
      )
    as.numeric(crps(pec)) # CRPS to extract IBS
    
  }, error = function(e) {
    NA_real_
    #   pec(
    #     object = list('coxph' = fit),
    #     formula = Surv(time, status) ~ 1,
    #     data = df,
    #     times = times,
    #     exact = F,
    #     reference = F,
    #     cens.model = "marginal",
    #     keep.matrix = F
    # )
  })
  return(ibs)
}

cox_metrics <- function(X, y, fit, k_times = 20) {
  cindex <- summary(fit)$concordance[1]
  ibs <- get_IBS_score(
    df = data.frame(cbind(X, y)),
    fit = fit,
    k_times = k_times
  )
  list(cindex = cindex, ibs = ibs)
}

## ---- Save results ---------------------------------
metric_log_init <- function() {
  tibble::tibble(
    batch_type = character(), # "BEXXAssoxx"
    norm_type = character(),  # "normXX"
    data_type = character(),
    n_train = integer(),
    iter = integer(),
    event_rate = numeric(),
    stratified = logical(),
    model = character(),     # "oracle-linear", "oracle-nl", "lasso"
    dataset = character(),   # "train" / "test"
    cind = numeric(),
    ibs = numeric(),
    lambda = numeric()       # NA unless lasso
  )
}

metric_log_add <- function(
  log,
  batchNormType, data_type, iter_i, train_size, ev_rate, stratify, model, dataset, metrics, lambda = NA_real_,
  out_dir = NULL # provide meta results file path if want to append per iteration
) {
  parts = stringr::str_split(batchNormType, "_", simplify = T)
  new_row <- tibble::tibble(
      batch_type=parts[1], 
      norm_type=parts[2],
      data_type=data_type,
      n_train = train_size,
      iter = iter_i,
      event_rate = ev_rate,
      stratified = stratify,
      model, dataset,
      cind = metrics$cindex, ibs = metrics$ibs, lambda
    )
  # append one line each function call
  if (!is.null(out_dir)) {
    dir.create(dirname(out_dir), showWarnings = FALSE, recursive = TRUE)
    write.table(
      new_row, file = out_dir, row.names = FALSE,
      col.names = !file.exists(out_dir),  # write header only once
      append = file.exists(out_dir),
      quote = FALSE, sep = "\t"
    )
  }
  dplyr::bind_rows(log, new_row)
}


# ---- Main function ----------------------------------------------------

# N = 20000; P = 538; n_batch = 15
# test_size = 1000
# max_train_size = 10000
# rho_cutoff = 0.9
# sim.dataType = stringr::str_split(dataType_list[1], "-p")[[1]][1]
# p = stringr::str_split(dataType_list[1], "-p")[[1]][2]
# train_sizes = c(100, 200, 500, 1000) #2000, 5000, 10000),
# splits_per_size = c(3, 5, 5, 5) #10, 10, 10),
# n_iter = 20
# he_train = 0
# he_test = 0
# beta_sort_train = 0
# beta_sort_test = 0
# norm_type = 0
# plot_km = T
# save_surv_data = T
# save_gene_data = T
# run_analysis = T
# stratify = 1
# save_results = F
# surv_folder = "sim_surv"
# results_file = "model_results_p30_v2.csv"
sim.survdata <- function(
  sim.data,
  N = 20000, P = 538, n_batch = 15,
  test_size = 1000,
  max_train_size = 10000,
  train_sizes = c(100, 200, 500, 1000, 2000, 5000, 10000),
  splits_per_size = c(3, 5, 5, 5, 10, 10, 10),
  n_iter = 20,
  sim.dataType = NULL, # "linear-moderate", "nl-quadratic"
  p = NULL, # 2, 5, 10, 30, 60
  he_train = 0,
  he_test = 0,
  beta_sort_train = 0,
  beta_sort_test = 0,
  norm_type = 0,
  rho_cutoff = 0.9,
  plot_km = F,
  save_surv_data = T,
  save_gene_data = T,
  run_analysis = T,
  stratify = 1,
  save_results = T,
  surv_folder = "sim_surv",
  results_file = NULL,
  date = NULL
) {
  # get fixed chunk len per batch per iteration (max needed)
  max_train_per_batch <- ceiling(max_train_size / n_batch)

  # Check normalization type ---
  nm <- norm_map[[as.character(norm_type)]]
  if (is.null(nm)) stop("Unsupported norm_type: ", norm_type)

  # Code scenario name ---
  # e.g., "BE11Asso00_normTC" ({batchType}_norm{normType})
  convert_num_to_indicator <- function(x) {
    case_when(x > 0 ~ 1, x == 0 ~ 0, x < 0 ~ -1)
  }
  sort_train = convert_num_to_indicator(beta_sort_train)
  sort_test = convert_num_to_indicator(beta_sort_test)
  batchNormType <- glue::glue(
    "BE{he_train}{he_test}Asso{sort_train}{sort_test}_norm{nm$name}"
  )
  
  # Initialize results df and output directory
  sim.dataType.p <- glue::glue("{sim.dataType}-p{p}")

  results <- metric_log_init()
  if (run_analysis) {
    if (!save_results) out_dir <- NULL
    else {
      results.dir <- file.path("results", "models", batchNormType, sim.dataType.p, "cox")
      dir.create(results.dir, showWarnings = F, recursive = T)

      if (is.null(results_file)) {
        out_dir <- file.path(results.dir, glue::glue("model_results_all_{date}.csv"))
      } else {
        out_dir <- file.path(results.dir, results_file)
      }
    }
  }

  # Start Loop ----------------------------------------------------
  
  message(glue::glue("Starting simulation for {batchNormType}-{sim.dataType} (p={p})"))

  for (train_size in train_sizes) {
# train_size=500
    n_splits <- splits_per_size[which(train_sizes == train_size)]

    for (iter_i in seq_len(n_iter)) {
# iter_i = 1
      message(paste0("Train size: ", train_size, " (Iteration ", iter_i, ")\n"))

      ## 1. Prepare data split ------------------------------------

      split <- extract_split_ids(train_size, test_size, iter_i, max_train_per_batch, n_batch, N)
      train.ids <- split[['train_idx']]
      test.ids <- split[['test_idx']]
      train_counts <- split[['train_counts']]
      test_counts <- split[['test_counts']]

      # Separate clean and batch-contaminated data
      ## [NOTE] latest attempt has clean data first, dirty data second!
      xtrue.train <- sim.data[train.ids, 1:P]
      xtrue.test <- sim.data[test.ids, 1:P]
      batch.train <- sim.data[train.ids, (P + 1):(2 * P)]
      batch.test <- sim.data[test.ids, (P + 1):(2 * P)]

      # Align marker ids in clean and BE data
      colnames(xtrue.train) <- colnames(xtrue.test) <- 
        colnames(batch.train) <- colnames(batch.test) <- g.names

      ## 2. Get simulated survival outcome ----------------------------

      ## [TEMPORARY SURV FOLDER!!]
      surv.file <- file.path("data", surv_folder, glue::glue("simSurv_{sim.dataType.p}_full.csv"))
      surv.out <- read.csv(surv.file)

      t.train = surv.out[as.character(train.ids), 'time']
      t.test = surv.out[as.character(test.ids), 'time']
      status.train = surv.out[as.character(train.ids), 'status']
      status.test = surv.out[as.character(test.ids), 'status']

      # 3. Insert batch effects ----------------------------------------------

      batch.id.train = rep(1:n_batch, times = train_counts)
      batch.id.test = rep(1:n_batch, times = test_counts)
      batch.id.unique = seq_len(n_batch)
      # reset index
      rownames(xtrue.train) = rownames(batch.train) = 1:train_size #names(t.train) = names(status.train) = 1:train_size
      rownames(xtrue.test) = rownames(batch.test) = 1:test_size #names(t.test) = names(status.test ) = 1:test_size

      ### HE Train ------------
      id = id.new = seq_len(train_size)

      if (he_train == 0) {
        x.train.count = xtrue.train
        x.train = log(x.train.count + 0.5)
      } else {
        ### NOTE: extract batch on the natural log scale
        xtrue.train.log = log(xtrue.train + 0.5)
        batch.train.log = log(batch.train + 0.5)
        ### [update 7/16] scale batch effects by 1/2
        ### [update 2/24/26] changed to normal scale
        batch.obs = (batch.train.log - xtrue.train.log)

        ### [Per Andy 7/21/2025]
        ### The batch effects should be sorted for all genes
        ### at the same time to preserve the correlation structure of the batch effects across genes.
        ### You can calculate the median of the median batch effects across genes in each batch
        ### and then simply “hard” sort the medians across batches. There is no need for partial sorting at this stage.

        ## get per batch median
        batch.obs$batch_id = batch.id.train
        batch.order = cbind(
          Median = apply(
            batch.obs %>%
              group_by(batch_id) %>%
              summarise(across(everything(), median)) %>%
              select(-batch_id),
            1, median),
          batch_id = seq(1, n_batch)
        ) %>%
          as.data.frame() %>%
          arrange(desc(Median)) %>%
          pull(batch_id)

        batch.obs$batch_id = factor(batch.obs$batch_id, levels = batch.order)
        batch.obs.sorted = batch.obs[order(batch.obs$batch_id), ] %>%
          select(-batch_id)

        ## sort batch by survival time
        if (beta_sort_train != 0) {
          temp = exp(beta_sort_train * log(t.train))
          id2 = seq_len(train_size)
          id.new = rep(NA, train_size)

          for (i in 1:(train_size - 1)) {
            probs = temp / sum(temp)
            sel = rmultinom(1, 1, probs)
            id.new[i] = id2[sel == 1]
            id2 = id2[sel != 1]
            temp = temp[sel != 1]
          }
          id.new[train_size] = id[!(id %in% id.new)]

          # Sort batch effect
          ### [Per Andy 7/21/2025
          ### Should be “batch.train.log = batch.obs.sorted+ xtrue.train.log[id.new,]”,
          ### i.e. it is the xtrue.train.log that needs to be partially sorted by id.new, but batch.obs.sorted.
          ### This is because in line 618 you sort t.train by id.new.
          batch.train.log = xtrue.train.log[id.new, ] + batch.obs.sorted
          batch.train = exp(as.data.frame(lapply(
            batch.train.log,
            cap_and_add_noise,
            cap = 19
          )))
        }

        # Final train data
        x.train.count = batch.train
        x.train = log(x.train.count + 0.5) # on natural log scale
      }

      ### HE Test -------------
      id.test = id.new.test = seq_len(test_size)

      if (he_test == 0) {
        x.test.count = xtrue.test
        x.test = log(x.test.count + 0.5)
      } else {
        xtrue.test.log = log(xtrue.test + 0.5)
        batch.test.log = log(batch.test + 0.5)
        batch.obs = (batch.test.log - xtrue.test.log)

        ## per batch median
        batch.obs$batch_id = batch.id.test
        batch.order = cbind(
          Median = apply(
            batch.obs %>%
              group_by(batch_id) %>%
              summarise(across(everything(), median)) %>%
              select(-batch_id),
            1, median),
          batch_id = seq(1, n_batch)
        ) %>%
          as.data.frame() %>%
          arrange(desc(Median)) %>%
          pull(batch_id)

        batch.obs$batch_id = factor(batch.obs$batch_id, levels = batch.order)
        batch.obs.sorted = batch.obs[order(batch.obs$batch_id), ] %>%
          select(-batch_id)

        ## Sort batch by survival time
        if (beta_sort_test != 0) {
          temp = exp(beta_sort_test * log(t.test))
          id2 = seq_len(test_size)
          id.new.test = rep(NA, test_size)

          for (i in 1:(test_size - 1)) {
            probs = temp / sum(temp)
            sel = rmultinom(1, 1, probs)
            id.new.test[i] = id2[sel == 1]
            id2 = id2[sel != 1]
            temp = temp[sel != 1]
          }
          id.new.test[test_size] = id.test[!(id.test %in% id.new.test)]

          # Sort batch effect
          batch.test.log = xtrue.test.log[id.new.test, ] + batch.obs.sorted
          batch.test = exp(as.data.frame(lapply(
            batch.test.log,
            cap_and_add_noise,
            cap = 19
          )))
        }

        # Final test data
        x.test.count = batch.test
        x.test = log(x.test.count + 0.5)
      }

      # 4. Normalization ---------------------------
      
      if (!nm$name %in% c('None', 'Quantile')) {
        train_fit <- nm$fun(raw = t(round(x.train.count)))
        test_fit <- nm$fun(raw = t(round(x.test.count)))

        sf_train <- train_fit$scalingFactor
        sf_test <- test_fit$scalingFactor

        x.train <- apply_sf_then_log(x.train.count, sf_train)
        x.test <- apply_sf_then_log(x.test.count, sf_test)
      } 
      if (nm$name == "Quantile") {
        # Quantile normalization
        x.train = t(preprocessCore::normalize.quantiles(t(x.train)))
        x_quantile = sort(x.train[1, ])

        x.test = t(apply(x.test, 1, function(u) {
          x_quantile[rank(u)]
        }))
      }
      colnames(x.train) <- colnames(x.test) <- g.names

      # 5. Save data ------------------------------

      ## ---- Final train and test data ----
      sim.train <- data.frame(
        batch_id = factor(batch.id.train),
        sample_idx = train.ids, # <-- save the selected indices used to generate this row
        time = t.train[id.new],
        status = status.train[id.new],
        # x.train,
        check.names = FALSE
      )
      sim.test <- data.frame(
        batch_id = factor(batch.id.test),
        sample_idx = test.ids, # <-- save test indices too
        time = t.test[id.new.test],
        status = status.test[id.new.test],
        # x.test,
        check.names = FALSE
      )
      tr_ev_rate <- sum(status.train) / train_size
      te_ev_rate <- sum(status.test) / test_size
      message(paste("Train event rate:", tr_ev_rate, " | ", "Test event rate:", te_ev_rate))

      ## ---- K-M plot of simulated train/test survival outcome ----
      if (plot_km) {
        plot_KM(df = rbind(
            sim.train %>% dplyr::mutate(type = 'Train'),
            sim.test %>% dplyr::mutate(type = 'Test')
        ))
      }

      ## ---- Save survival + batch id + sample id  ----
      if (save_surv_data) {
        # output directory        
        surv_dir <- file.path(
          "data", batchNormType, sim.dataType.p, paste0("train", train_size))
        dir.create(surv_dir, showWarnings = FALSE, recursive = TRUE) 

        # consistent file stem for this experiment
        file_stem <- glue::glue(
          "{batchNormType}_{sim.dataType.p}_train{train_size}_iter{iter_i}"
        ) 
        write.csv(
          sim.train |> dplyr::select(-batch_id),
          file.path(surv_dir, glue::glue("simSurv_train_{file_stem}.csv")),
          row.names = FALSE
        )
        write.csv(
          sim.test |> dplyr::select(-batch_id),
          file.path(surv_dir, glue::glue("simSurv_test_{file_stem}.csv")),
          row.names = FALSE
        )
      }
      ## ---- Save partitioned + normalized gene expression ----
      if (save_gene_data) {
        # output director
        gene_tr_dir = file.path('data', "gene_exp", batchNormType, glue::glue("train{train_size}"))
        gene_te_dir = file.path('data', "gene_exp", batchNormType)
        dir.create(gene_tr_dir, showWarnings = F, recursive = T)

        x_train_file = glue::glue("simGeneExp_train_{batchNormType}_train{train_size}_iter{iter_i}.parquet")
        x_test_file = glue::glue("simGeneExp_test_{batchNormType}_iter{iter_i}.parquet")
        x_train_file_path = file.path(gene_tr_dir, x_train_file)
        x_test_file_path = file.path(gene_te_dir, x_test_file)

        if (!file.exists(x_train_file_path)) {
          write_parquet(
            as.data.frame(cbind(x.train, batch_id = batch.id.train)),
            x_train_file_path)
        }
        if (!file.exists(x_test_file_path)) {
          write_parquet(
            as.data.frame(cbind(x.test, batch_id = batch.id.test)), 
            x_test_file_path)
        }
      }

      ## 6. Analysis --------------------------------------------------
      
      if (run_analysis) { 

        ## * ---- Pre-screen -------------------------------------
        means_train=colMeans(x.train)
        xcorr_train=cor(x.train)
        exclude_train=c()
        for(i in 1:ncol(x.train)){
          for(j in i:ncol(x.train)){
            if (is.na(xcorr_train[i,j])) next
            if(xcorr_train[i,j] > rho_cutoff & xcorr_train[i,j]!=1){
              exclude_train_1=ifelse(means_train[g.names==rownames(xcorr_train)[i]]>=means_train[g.names==colnames(xcorr_train)[j]],
                            colnames(xcorr_train)[j], rownames(xcorr_train)[i])
              exclude_train=c(exclude_train, exclude_train_1)
            }
          }
        }
        # [ADD 3/5/26] remove 0 variance markers (corr=NA)
        sd0 <- names(which(apply(x.train, 2, sd, na.rm=TRUE) == 0))
        if (length(sd0)!=0) {exclude_train=c(exclude_train, sd0)}
        if(length(unique(exclude_train))==0){
          x_train=x.train
        }else{
          x_train=x.train[,!(colnames(x.train) %in% unique(exclude_train))]
        }

        x_test = x.test

        # select true markers for oracle analysis (read from presaved beta values)
        beta.file <- file.path("data", surv_folder, glue::glue("beta0-p{p}.csv"))
        selected.geneid <- read.csv(beta.file)[,1]

        exclude_true <- intersect(unique(exclude_train), selected.geneid)

        message(paste("\texcluded", length(unique(exclude_train)), "markers with corr >", rho_cutoff))
        # message(paste(unique(exclude_train), collapse = ", "))
        message(paste("\texcluded",length(exclude_true),"true markers:\n", paste(exclude_true, collapse = ", ")))

        x0_train = x.train[, selected.geneid]
        x0_test = x.test[, selected.geneid]
        y_train = sim.train[,c('time', 'status','batch_id')]
        y_test = sim.test[,c('time', 'status', 'batch_id')]

        ## * ---- Oracle linear (non-strat) ----------------------
        formula <- as.formula(glue::glue(
            "Surv(time, status) ~ {paste(colnames(x0_train), collapse = '+')}"
          ))
        
        coxfit_o <- fit_coxph_safe(formula, x0_train, y_train)
        metrics_o_tr <- cox_metrics(x0_train, y_train, coxfit_o)
        coxfit_o_test <- fit_coxph_safe(formula, x0_test, y_test, beta_init = summary(coxfit_o)$coefficient[, 1])
        metrics_o_te <- cox_metrics(x0_test, y_test, coxfit_o_test)
        
        ## * ---- Oracle nonlinear (optional) --------------------
        metrics_nl_tr <- metrics_o_tr
        metrics_nl_te <- metrics_o_te

        if (!grepl("linear", sim.dataType)) { #%in% c("linear-moderate", "linear-weak")) {

          if (sim.dataType=='nl-interaction') interact.gene <- "hsa.miR.33a.1"
          x0_t_train <- x0_transform(sim.dataType, x0_train, interact.gene)
          x0_t_test <- x0_transform(sim.dataType, x0_test, interact.gene)

          coxph_nl <- fit_coxph_safe(formula, x0_t_train, y_train)
          metrics_nl_tr <- cox_metrics(x0_t_train, y_train, coxph_nl)

          beta_nl <- summary(coxph_nl)$coefficient[, 1]
          coxph_nl_test <- fit_coxph_safe(formula, x0_t_test, y_test, beta_init = beta_nl)
          metrics_nl_te <- cox_metrics(x0_t_test,  y_test, coxph_nl_test)
        }
        
        ## * ---- Lasso (non-strat) ----------------------------------
        cv_l = cv.glmnet(
          as.matrix(x_train), # as.matrix(x_train[,lkhd_u >= lkhd_p_thres]),
          Surv(y_train$time, y_train$status),
          nfolds = n_splits,
          family = "cox",
          alpha = 1, # Lasso
          standardize = F
        )
        lambda_l = cv_l$lambda.min

        # Fit final lasso regression
        glmnet_l = glmnet(
          as.matrix(x_train), # as.matrix(x_train[,lkhd_u >= lkhd_p_thres]),
          Surv(y_train$time, y_train$status),
          family = "cox",
          alpha = 1,
          lambda = lambda_l,
          standardize = F
        )
        b_l_temp = as.vector(coef(glmnet_l))
        b_l_sel = b_l_temp[b_l_temp != 0]

        metrics_l_tr <- metrics_l_te <-list(cindex = NA_real_, ibs = NA_real_)   
        if (length(b_l_sel) > 0) {
          x_train_sel <- data.frame(x_train[, b_l_temp != 0, drop = FALSE])
          x_test_sel <- data.frame(x_test[, colnames(x_train_sel), drop = FALSE])
          formula_l <- as.formula(glue::glue(
            "Surv(time, status) ~ {paste(colnames(x_train_sel), collapse='+')}"
          ))
          coxph_l <- fit_coxph_safe(formula_l, x_train_sel, y_train, beta_init = b_l_sel)
          metrics_l_tr <- cox_metrics( x_train_sel, y_train,coxph_l)
          coxph_l_test <- fit_coxph_safe(formula_l, x_test_sel, y_test, beta_init = b_l_sel)
          metrics_l_te <- cox_metrics( x_test_sel, y_test,coxph_l_test)
        }

        if (stratify) {

          ## * ---- Oracle linear (strat) --------------------------
          formula <- as.formula(paste(
            "Surv(time, status) ~", paste(colnames(x0_train), collapse = "+"), "+ strata(batch_id)"
          ))
          coxph_os <- fit_coxph_safe(formula, x0_train, y_train, stratify = T)
          metrics_os_tr <- cox_metrics(x0_train, y_train, coxph_os)
          coxph_os_test <- fit_coxph_safe(
            formula, x0_test, y_test, 
            beta_init=summary(coxph_os)$coefficient[, 1], 
            stratify = T)
          metrics_os_te <- cox_metrics(x0_test, y_test, coxph_os_test)

          ## * ---- Oracle nonlinear (strat) -------------------------
          metrics_nls_tr <- metrics_os_tr
          metrics_nls_te <- metrics_os_te #list(cindex = NA_real_, ibs = NA_real_)
          
          if (!grepl("linear", sim.dataType)) { #%in% c("linear-moderate", "linear-weak")) {
            x0_t_train <- x0_transform(sim.dataType, x0_train, interact.gene = interact.gene)
            x0_t_test <- x0_transform(sim.dataType, x0_train, interact.gene = interact.gene)

            coxph_nls <- fit_coxph_safe(formula, x0_t_train, y_train, stratify = stratify)
            metrics_nls_tr <- cox_metrics(x0_t_train, y_train, coxph_nls)

            beta_nls <- summary(coxph_nls)$coefficient[, 1]
            coxph_nls_test <- fit_coxph_safe(formula, x0_t_test, y_test, beta_init = beta_nls, stratify = stratify)
            metrics_nls_te <- cox_metrics(x0_t_test, y_test, coxph_nl_test)
          }

          ## * ---- LASSO (strat with manual CV) ----------------------
          folds <- caret::createMultiFolds(
            y_train$batch_id,
            k = n_splits,
            times = 1
          )
          cv_cindex = c()
          lasso_final = NULL
          for (j in seq_along(folds)) {
            fold_idx = folds[[j]]
            y_fold_train = y_train[fold_idx, ]
            y_fold_valid = y_train[-fold_idx, ]
            x_fold_train = x_train[fold_idx, ]
            x_fold_valid = x_train[-fold_idx, ]

            val_cindex <- tryCatch({# Attempt to fit glmnet Cox model
              lasso = glmnet(
                as.matrix(x_fold_train),
                stratifySurv(
                  Surv(y_fold_train$time, y_fold_train$status),
                  strata = y_fold_train$batch_id),
                family = 'cox',
                alpha = 1,
                nlambda = 20,
                standardize = F
              )
              if (is.null(lasso_final)) {lasso_final = lasso}
              val_preds = predict(lasso, newx = as.matrix(x_fold_valid), type = 'link')
              apply(val_preds, 2, function(x) {
                survcomp::concordance.index(x, y_fold_valid$time, y_fold_valid$status)$c.index
              })
            },
            error = function(e) {
              message(sprintf("⚠Fold %d failed: %s", j, e$message))
              rep(NA, 20) # same length as nlambda
            })
            cv_cindex = cbind(cv_cindex, val_cindex)
          }

          metrics_ls_tr <- metrics_ls_te <- list(cindex = NA_real_, ibs = NA_real_)   
          if (!is.null(lasso_final) && any(colSums(!is.na(cv_cindex)) > 0)) {
            lambda_ls = lasso_final$lambda[which.max(rowMeans(cv_cindex, na.rm = T))]

            glmnet_ls = glmnet(
              as.matrix(x_train), # as.matrix(x_train[,lkhd_u >= lkhd_p_thres]),
              stratifySurv(
                Surv(y_train$time, y_train$status),
                y_train$batch_id),
              family = "cox",
              alpha = 1,
              lambda = lambda_ls,
              standardize = F
            )
            b_ls_temp = as.vector(coef(glmnet_ls))
            b_ls_sel = b_ls_temp[b_ls_temp != 0]

            if (length(b_ls_sel) > 0) {
              x_train_sel = data.frame(x_train[, b_ls_temp != 0, drop = FALSE])
              x_test_sel = data.frame(x_test[, colnames(x_train_sel), drop = FALSE])
              formula_ls <- as.formula(glue::glue(
                "Surv(time, status) ~ {paste(colnames(x_train_sel), collapse='+')} + strata(batch_id)"
              ))
              coxph_ls <- fit_coxph_safe(formula_ls, x_train_sel, y_train, beta_init = b_ls_sel)
              metrics_ls_tr <- cox_metrics(x_train_sel, y_train, coxph_ls)
              coxph_ls_test <- fit_coxph_safe(formula_ls, x_test_sel, y_test, beta_init = b_ls_sel)
              metrics_ls_te <- cox_metrics(x_test_sel, y_test,coxph_ls_test)
            }
          } else {
            message("⚠All folds failed or no valid C-index values — setting lambda_min = NA")
            lambda_ls = NA
          }
        } # end of stratified analysis

        ## * ---- Save model results ------------
        # Log results per iteration
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, tr_ev_rate, stratify=F, "oracle-linear", "train", metrics_o_tr, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, te_ev_rate, stratify=F, "oracle-linear", "test", metrics_o_te, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, tr_ev_rate, stratify=F, "oracle-nl", "train", metrics_nl_tr, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, te_ev_rate, stratify=F, "oracle-nl", "test", metrics_nl_te, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, tr_ev_rate, stratify=F, "lasso", "train", metrics_l_tr, lambda_l, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, te_ev_rate, stratify=F, "lasso", "test", metrics_l_te, lambda_l, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, tr_ev_rate, stratify=T, "oracle-linear", "train", metrics_os_tr, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, te_ev_rate, stratify=T, "oracle-linear", "test", metrics_os_te, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, tr_ev_rate, stratify=T, "oracle-nl", "train", metrics_nls_tr, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, te_ev_rate, stratify=T, "oracle-nl", "test", metrics_nls_te, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, tr_ev_rate, stratify=T, "lasso", "train", metrics_ls_tr, lambda_ls, out_dir=out_dir)
        results <- metric_log_add(results, batchNormType, sim.dataType.p, iter_i, train_size, te_ev_rate, stratify=T, "lasso", "test", metrics_ls_te, lambda_ls, out_dir=out_dir)

      } # end of run_analysis
    } # end of all iter for loop
  } # end of all subsets for loop

  return(results)
}


