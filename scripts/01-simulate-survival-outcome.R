#################################################
# Date      Note
# 17MAR2026 Generate surv data with 2x beta coefs
# TODO - check var(h.log) to resolve weak signal being more predictive
# The variance of h.log ~ sum(beta_j^2 * var_j). 
# So moderate has more variance in the linear predictor->larger spread in survival risk
# But we'd expect moderate to be more discriminative
# otherwise harder to pick up signal for prediction.
#################################################
library(nanoparquet)
library(ggsurvfit)
library(edgeR)


# Load Raw Augmented Gene Counts ---------------------------------------------
P = 538
g.names = read.csv("data/augmentation_miRNA_names_538.csv")$gene # gene names
nonzero_position = c(56, 57, 79, 100, 112, 141, 184, 295, 466, 518) # RWD true gene positions
#  [1] "hsa.miR.1277.3p.1" "hsa.miR.1277.5p.1" "hsa.miR.133a.2..1" "hsa.miR.144..1"    "hsa.miR.148b..1"   
#  [6] "hsa.miR.181c..1"   "hsa.miR.204.1"     "hsa.miR.33a.1"     "hsa.miR.598.1"     "hsa.miR.887.1"  

nonzero.genes.list = list(
  '2'  = c("hsa.miR.181c..1", "hsa.miR.33a.1"),
  '5'  = c("hsa.miR.181c..1", "hsa.miR.33a.1", "hsa.miR.204.1", "hsa.miR.598.1", "hsa.miR.887.1"),
  '10' = g.names[nonzero_position],
  '30' = c(g.names[nonzero_position],
          "hsa.miR.183.1", "hsa.miR.129.2.3p.1", "hsa.miR.509.3p.3..1", "hsa.miR.370..1", "hsa.miR.182.1",      
          "hsa.miR.381..1", "hsa.miR.134.1", "hsa.miR.105.2..1", "hsa.miR.196b..1", "hsa.miR.34b.1",
          "hsa.miR.767.1", "hsa.miR.199b.5p.1", "hsa.miR.99a..1", "hsa.miR.675.1", "hsa.miR.499.1",
          "hsa.miR.876.5p.1", "hsa.miR.488.1", "hsa.miR.3117.1", "hsa.miR.196b.1", "hsa.miR.758..1")
  # ----- old -----
  # '30' = c(g.names[nonzero_position],
  #         'hsa.miR.493.1','hsa.miR.582.1', 'hsa.miR.299.5p.1', 'hsa.miR.374a..1','hsa.miR.889.1',
  #         'hsa.miR.410.1', 'hsa.miR.382.3p.1', 'hsa.miR.494.1', 'hsa.miR.335.3p.1', 'hsa.miR.409.5p.1',
  #         'hsa.miR.421.1','hsa.miR.337.1', 'hsa.miR.150.1', 'hsa.miR.382.5p.1', 'hsa.miR.340.1',
  #         'hsa.miR.495.1','hsa.miR.22..1', 'hsa.miR.299.3p.1', 'hsa.miR.532.3p.1'),
  # '60' = c("hsa.miR.1277.3p.1", "hsa.miR.1277.5p.1", "hsa.miR.133a.2..1", "hsa.miR.144..1",
  #         "hsa.miR.148b..1", "hsa.miR.181c..1", "hsa.miR.204.1", "hsa.miR.33a.1",
  #         "hsa.miR.598.1", "hsa.miR.887.1", "hsa.miR.493.1", "hsa.miR.582.1",
  #         "hsa.miR.299.5p.1", "hsa.miR.374a..1", "hsa.miR.889.1", "hsa.miR.410.1",
  #         "hsa.miR.382.3p.1", "hsa.miR.33b.1", "hsa.miR.494.1", "hsa.miR.335.3p.1",
  #         "hsa.miR.409.5p.1", "hsa.miR.421.1", "hsa.miR.337.1", "hsa.miR.150.1",
  #         "hsa.miR.382.5p.1", "hsa.miR.340.1", "hsa.miR.495.1", "hsa.miR.22..1",
  #         "hsa.miR.299.3p.1", "hsa.miR.532.3p.1", "hsa.miR.483.5p.1", "hsa.miR.31.3p.1",
  #         "hsa.miR.376a.1.5p.1", "hsa.miR.654.1", "hsa.miR.188.1", "hsa.miR.154.3p.1",
  #         "hsa.miR.654..1", "hsa.miR.95.1", "hsa.miR.148a..1", "hsa.miR.501.3p.1",
  #         "hsa.miR.758.1", "hsa.miR.455.5p.1", "hsa.miR.181d.1", "hsa.miR.301b.1",
  #         "hsa.miR.221..1", "hsa.miR.362.3p.1", "hsa.miR.493..1", "hsa.miR.432.1",
  #         "hsa.miR.708.1", "hsa.miR.483.3p.1", "hsa.miR.370.1", "hsa.miR.183.1",
  #         "hsa.miR.487a.3p.1", "hsa.miR.133b.1", "hsa.miR.543.1", "hsa.miR.656.1",
  #         "hsa.miR.1185.5p.2..1", "hsa.miR.499.1", "hsa.miR.99a..1", "hsa.miR.196a.2..2")         
)


parquet_file_path <- file.path("data", "simGeneExp_full_preprocessed.parquet")

if (! file.exists(parquet_file_path)) {
  # load raw augmented data (attempt 18FEB2026)
  sim.data = read.csv(
    file.path("raw-data",
      "with_without_15_batch_CVAE_20000perbatch_02182026.csv"
    ), header = T)[, 1:(P * 2)]

  ## Cap extremely large values
  ## [update 6/12/2025] Cap each marker (on the log2 scale) by (max + 0.1*sd)
  ## estimated from the 108 real-world sarcoma samples; then add Gaussian noise
  ## to the log2-transformed counts with mean=0 and standard dev = 0.01*sd
  cap_and_add_noise <- function(x, cap, sd) {
    x[x < 0] = 0 # Check if there's any negative values
    new_cap = cap + 0.1 * sd
    tryCatch({
      xs = x[is.infinite(x) | (x > new_cap)]
      xs_hat = pmin(xs, new_cap + rnorm(length(xs), 0, 0.01 * sd))
      x[is.infinite(x) | (x > new_cap)] = xs_hat
      return(x)
    }, 
    error = function(e) {
      return(x)
    })
  }
  ## [Note] DESeq normalization has a built-in as.integer() processing that takes
  ## maximum of ~exp(21) or else error would occure. Therefore we lower the
  ## cap value to 19 (empirically chosen) when normalization method is DEseq
  ## [NOTE] newest augmented attempt has no extremely large values but kept this step for rigor
  sim.data <- exp(as.data.frame(lapply(
    log(sim.data + 1),
    cap_and_add_noise,
    cap = 19)
  )) - 1
  # Save capped dataset (parquet for fast loading)
  write_parquet(
    sim.data, 
    file.path("data", "simGeneExp_full_preprocessed.parquet")
  )
} else {
  sim.data = read_parquet(file.path("data", "simGeneExp_full_preprocessed.parquet"))
}


## [NOTE] attempt 18FEB2026: clean data first, dirty data last 
x.true <- sim.data[, 1:P]
colnames(x.true) <- g.names


# # * ---------- FOR TEST PURPOSES ONLY ! ---------- * 
# x.true.log <- log(x.true + 0.5)
# vars_2 <- sapply(x.true.log[,nonzero.genes.list[['2']]], var)
# vars_5 <- sapply(x.true.log[,nonzero.genes.list[['5']]], var)
# vars_10 <- sapply(x.true.log[,nonzero.genes.list[['10']]], var)
# vars_30 <- sapply(x.true.log[,nonzero.genes.list[['30']][! nonzero.genes.list[['30']] %in% nonzero.genes.list[['10']]]], var)
# # MANUALLY SELECT another 20 markers with 1.98 variance
# # all_vars <- sapply(x.true.log[,g.names[!g.names %in% nonzero.genes.list[['30']]]], var)  
# # sel_vars_20 <- sample(all_vars[all_vars >= 1.8 & all_vars <= 2.2], 20, replace=F)
# # sel_genes_20 <- names(sel_vars_20)
# vars_30_2 <- sapply(x.true.log[,nonzero.genes.list[['30_2']][! nonzero.genes.list[['30_2']] %in% nonzero.genes.list[['10']]]], var)

# cat("mean (p=2):", mean(vars_2), "\n")
# cat("mean (p=5):", mean(vars_5), "\n")
# cat("mean (p=10):", mean(vars_10), "\n")
# cat("mean (p=30; excluding 10 true markers):", mean(vars_30), "\n")
# cat("mean (p=30; v2; excluding 10 true markers):", mean(vars_30_2), "\n")

# surv_path = file.path('data', 'sim_surv')
# beta_2 <-  read.csv(list.files(surv_path, pattern = glue::glue("beta0-p2.csv"), full.names = T))
# beta_5 <-  read.csv(list.files(surv_path, pattern = glue::glue("beta0-p5.csv"), full.names = T))
# beta_10 <- read.csv(list.files(surv_path, pattern = glue::glue("beta0-p10.csv"), full.names = T))
# beta_30 <- read.csv(list.files(surv_path, pattern = glue::glue("beta0-p30.csv"), full.names = T))
# beta_30_2 <- read.csv(list.files(surv_path, pattern = glue::glue("beta0-p30_2.csv"), full.names = T))

# cat("sum abs(beta) (p=2):", sum(abs(beta_2$beta0)), "\n")
# cat("sum abs(beta) (p=3):", sum(abs(beta_5$beta0)), "\n")
# cat("sum abs(beta) (p=10):", sum(abs(beta_10$beta0)), "\n")
# cat("sum abs(beta) (p=30; excluding 10 true markers):", sum(abs(beta_30$beta0)), "\n")
# cat("sum abs(beta) (p=30; v2; excluding 10 true markers):", sum(abs(beta_30_2$beta0)), "\n")
# * ------------------------------------------------ *


# Simulate survival outcome ------------------------------------------------

simulate_T <- function(h0, h.log, n, max.censor.time) {
  # simulate uncensored survival time
  t0 = -log(runif(n)) * h0 / exp(h.log)
  # create censoring time
  ct = runif(n, 0, max.censor.time) ## capped by a max time
  t = ifelse(t0 <= ct, t0, ct)
  delta = 1 * (t0 <= ct)
  ncase = sum(delta)
  event.rate = ncase/ n
  cat(glue::glue("Event rate: {event.rate}\n\n")) ## about 75% event rate
  # print(summary(t)) ## summarize simulated survival time

  return(list(time = t, status = delta, ev_rate = event.rate))
}


calculate_beta0 <- function(x.log, p, scale = 1) {
  sd.genes = apply(x.log, 1, sd)
  beta_scale = 10 / p * scale
  rep(c(1, -1), ceiling(p / 2))[1:p] * round(beta_scale / sd.genes, 3)
}


# sim.dataType = "linear"
# max.censor.time = 200
# p = "30_2"
# nonzero.genes = nonzero.genes.list[[as.character(p)]]
# h0 = 1.5
# save_surv_data = F
# surv_folder = "sim_surv"
# plot_km = T
# beta.scale = NULL
simulate_survival <- function(
  x.true,
  sim.dataType, 
  h0, 
  p = 10,
  beta.scale = NULL,
  nonzero.genes = NULL,
  max.censor.time = 200,
  save_surv_data = T,
  save_beta = F,
  surv_folder = "sim_surv",
  plot_km = F) {

  set.seed(1234)

  surv_path <- file.path("data", surv_folder)
  if (is.null(beta.scale)) {
    beta.scale <- 1
    suffix <- glue::glue("p{p}")
  } else {
    suffix <- glue::glue("p{p}-{beta.scale}x")
  }
  beta_file <- list.files(surv_path, pattern = glue::glue("beta0-{suffix}.csv"), full.names = T)

  if (length(beta_file) != 0) {
    beta.geneid <- read.csv(beta_file)
    selected.geneid <- beta.geneid$X
    beta0 <- beta.geneid$beta0
    workingt <- t(x.true)

  } else {
    # ----------- Select true markers from data subset ------------
    workingt = t(x.true)
    ## Filter genes by count per million
    keep = (rowSums(cpm(workingt) > 2) >= 10) # adjustable
    workingt = workingt[keep, ]
    ## Select genes with higher TC and variance
    tc = rowSums(workingt)
    mid.tc = (tc > median(tc) & tc < quantile(tc, 0.75)) # select genes with total count in the upper quartile
    cv = apply(workingt, 1, sd) / rowMeans(workingt)
    mid.cv = (cv > quantile(cv, 0.5)) # select genes with higher variance (SD/mean)
    candidate.geneid = rownames(workingt)[mid.tc & mid.cv]
    
    if (is.null(nonzero.genes)) {
      # Select p (default=10) true genes associated with survival
      if (p <= length(candidate.geneid)) {
        p.vec = names(nonzero.genes.list)[names(nonzero.genes.list) < p]
        p.sub = p.vec[length(p.vec)]
        geneid.sub = nonzero.genes.list[[as.character(p.sub)]]
        geneid.add = sample(
          candidate.geneid[!candidate.geneid %in% geneid.sub], p - p.sub, replace = F
        )
        selected.geneid = c(geneid.sub, geneid.add)
        # selected.geneid = sample(candidate.geneid, p, replace = F)
      } else {
        selected.geneid = candidate.geneid
      }
    } else {
      selected.geneid = nonzero.genes
    }
    beta0 <- NULL
  }

  # ------ Simulate survival time with transformed marker data ------

  # log-transform selected true genes
  xtrue <- workingt[selected.geneid, ]
  xtrue.log <- log(xtrue + 0.5)
  n_total <- ncol(xtrue)

  ## linear signals --------
  if (sim.dataType == 'linear') {
    # sim.dataType = glue::glue('linear-p{p}') #stringr::str_split(sim.dataType, '-', simplify = T)[2]
    cat(glue::glue("Simulating from linear function with {p} true markers\n"))

    if (is.null(beta0)) beta0 <- calculate_beta0(xtrue.log, p, scale = beta.scale)
    
    ## calculate log risk
    h.log = t(xtrue.log) %*% beta0

  ## 10 quadratic terms -----------------
  } else if (sim.dataType == 'nl-quadratic') {
    cat("Simulating with quadratic transformation\n")

    ## Scale marker expression data
    xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
    xtrue.nl.log = as.matrix(xtrue.log.scaled)^2

    if (is.null(beta0)) beta0 <- calculate_beta0(xtrue.nl.log, p, scale = beta.scale)
    h.log = t(xtrue.nl.log) %*% beta0
    
  ## 10 interaction terms -----------------
  } else if (sim.dataType == 'nl-interaction') {
    cat("Simulating with interaction effects\n")

    interact.gene = sample(selected.geneid, 1)
    cat(glue::glue("Interaction term: {interact.gene}\n\n"))

    xtrue.log.scaled = t(apply(xtrue.log, 1, scale))
    interact.gene.val = xtrue.log.scaled[interact.gene, ]
    xtrue.nl.log = sweep(xtrue.log.scaled, 2, interact.gene.val, FUN = "*")

    if (is.null(beta0)) beta0 <- calculate_beta0(xtrue.nl.log, p, scale = beta.scale)
    h.log = t(xtrue.nl.log) %*% beta0

  ## 10 sine terms ----------------------
  } else if (sim.dataType == 'nl-sine') {
    cat("Simulating with non-linear sine-transformed terms\n\n")
    xtrue.nl.log = t(sin(t(xtrue.log))) # sine transformation
    if (is.null(beta0)) beta0 <- calculate_beta0(xtrue.nl.log, p, scale = beta.scale)
    h.log = t(xtrue.nl.log) %*% beta0
  }

  hist(h.log)
  message(paste("Variance of log hazard:", var(h.log)))
  message('Selected true markers:\n', paste(selected.geneid, collapse = ", "))
  message('True coefficients:\n', paste(beta0, collapse = ", "))

  ## Simulated survival time and censoring outcome
  surv.obj = simulate_T(h0, h.log, n_total, max.censor.time) 
  ev_rate = surv.obj$ev_rate
  surv.out = data.frame(time=surv.obj$time, status=surv.obj$status)
  rownames(surv.out) <- colnames(workingt)

  ## ---- 1) Save survival + batch id + sample id  ----
  if (save_surv_data) {
    dir.create(surv_path, recursive = T, showWarnings = F)
    write.csv(surv.out, 
      file.path(surv_path, glue::glue("simSurv_{sim.dataType}-{suffix}_full.csv")),
      row.names = F
    )
  }
  ## ---- 2) Save true beta coefs ----
  if (save_beta) {
    dir.create(surv_path, recursive = T, showWarnings = F)
    write.csv(data.frame(beta0),
      file.path(surv_path, glue::glue("beta0-{suffix}.csv"))
    )
  }
  out <- list(
    surv = surv.out,
    event_rate = ev_rate,
    p = NULL
  )
  if (plot_km) {
    p.km <- survfit(
        Surv(time, status) ~ 1,
        data = surv.out,
        conf.type = "log-log"
      ) |>
      ggsurvfit() +
      labs(title = glue::glue('Survival for {sim.dataType} association (p={p})'), x = 'Time', y = 'Survival Probability') +
      coord_cartesian(xlim = c(0, max.censor.time), ylim = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, max.censor.time, by = 24)) +
      add_confidence_interval() +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      theme_bw()
    print(p.km)
    out$p <- p.km
  }
  out
}


## linear with varying signal strength ---------------------------------------------------------------------------------------
# out = simulate_survival(x.true, 'linear', h0=2e-2, p=2, nonzero.genes=nonzero.genes.list[['2']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'linear', h0=5e5,   p=5, beta.scale=2, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'linear', h0=5.9e3, p=5, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'linear', h0=1.5e-2,p=10, nonzero.genes=nonzero.genes.list[['10']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'linear', h0=9.5,   p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
# out = simulate_survival(x.true, 'linear', h0=8e-3, p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = F, surv_folder = "sim_surv")
# out = simulate_survival(x.true, 'linear', h0=0.3, p=60, nonzero.genes=nonzero.genes.list[['60']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")


## Nonlinear (quadratic) with varying signal strength ---------------------------------------------------------------------------------------
# out = simulate_survival(x.true, 'nl-quadratic', h0=10, p=2, nonzero.genes=nonzero.genes.list[['2']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-quadratic', h0=37.5, p=5, beta.scale=2, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-quadratic', h0=55.5, p=5, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-quadratic', h0=8,    p=10, nonzero.genes=nonzero.genes.list[['10']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-quadratic', h0=20.5, p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
# out = simulate_survival(x.true, 'nl-quadratic', h0=9.75, p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
# out = simulate_survival(x.true, 'nl-quadratic', h0=18.5, p=60, nonzero.genes=nonzero.genes.list[['60']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")


## Nonlinear (interaction) with varying signal strength ---------------------------------------------------------------------------------------
# out = simulate_survival(x.true, 'nl-interaction', h0=7, p=2, nonzero.genes=nonzero.genes.list[['2']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-interaction', h0=39, p=5, beta.scale=2, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-interaction', h0=48, p=5, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-interaction', h0=16, p=10, nonzero.genes=nonzero.genes.list[['10']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-interaction', h0=40, p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
# out = simulate_survival(x.true, 'nl-interaction', h0=15, p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")


## Nonlinear (sine) with varying signal strength ---------------------------------------------------------------------------------------
# out = simulate_survival(x.true, 'nl-sine', h0=2, p=2, nonzero.genes=nonzero.genes.list[['2']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-sine', h0=.85, p=5, beta.scale=2, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-sine', h0=8.5, p=5, nonzero.genes=nonzero.genes.list[['5']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-sine', h0=115, p=10, nonzero.genes=nonzero.genes.list[['10']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
out = simulate_survival(x.true, 'nl-sine', h0=37,  p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")
# out = simulate_survival(x.true, 'nl-sine', h0=125, p=30, nonzero.genes=nonzero.genes.list[['30']], plot_km = T, save_surv_data = T, surv_folder = "sim_surv")


# ----------- * Test purpose * -----------
## P = 10, Beta scale = 1 --------------------------------------------------------------------------------------
# out = simulate_survival(x.true, 'linear-moderate', h0=, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'linear-weak', h0=, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'nl-quadratic', h0=, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'nl-interaction', h0=, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'nl-sine', h0=, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")

## P = 10, Beta scale = 2 --------------------------------------------------------------------------------------
# out = simulate_survival(x.true, 'linear-moderate', h0=3e-6, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'linear-weak', h0=5.2e-7, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'nl-quadratic', h0=10, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'nl-interaction', h0=21, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")
# out = simulate_survival(x.true, 'nl-sine', h0=155, beta_scale = beta_scale, nonzero.genes = nonzero.genes, plot_km = T, save_surv_data = F, surv_folder = "sim_surv_2")

## P = 2, Beta scale = 5 ---------------------------------------------------------------------------------------
# out = simulate_survival(x.true, 'linear-moderate', h0=2e-2, p=2, beta_scale=5, nonzero.genes=nonzero.genes, plot_km = T, save_surv_data = T, surv_folder = "sim_surv_p_2")
# out = simulate_survival(x.true, 'linear-weak', h0=6.25, p=2, beta_scale=5, nonzero.genes=nonzero.genes, plot_km = T, save_surv_data = T, surv_folder = "sim_surv_p_2")
