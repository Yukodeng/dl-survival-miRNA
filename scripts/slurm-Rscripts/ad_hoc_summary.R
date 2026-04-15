summarize_block <- function(results, n_train) {
  get_mean <- function(strat, mdl, ds) {
    results %>%
      filter(n_train == !!n_train, stratified == strat, model == mdl, dataset == ds) %>%
      summarise(c_m = mean(cind, na.rm=T), ibs_m = mean(ibs, na.rm=T)) %>%
      pull(c_m) %>%
      { ifelse(is.nan(.), NA_real_, .) }
  }
  fmt <- function(x) if (is.na(x)) "NA" else formatC(x, digits = 3, format = "f")

  # Non-stratified
  o_tr  <- fmt(get_mean(FALSE, "oracle-linear", "train")); o_te  <- fmt(get_mean(FALSE, "oracle-linear", "test"))
  nl_tr <- fmt(get_mean(FALSE, "oracle-nl",     "train")); nl_te <- fmt(get_mean(FALSE, "oracle-nl",     "test"))
  l_tr  <- fmt(get_mean(FALSE, "lasso",         "train")); l_te  <- fmt(get_mean(FALSE, "lasso",         "test"))

  out <- glue::glue(
"=============== Non-stratified ===============
(Oracle linear)     Train: {o_tr} |  Test: {o_te}
(Oracle nonlinear)  Train: {nl_tr} |  Test: {nl_te}
(Lasso)             Train: {l_tr} |  Test: {l_te}"
  )
  if (any(results$stratified) ) {
    os_tr  <- fmt(get_mean(TRUE, "oracle-linear", "train")); os_te  <- fmt(get_mean(TRUE, "oracle-linear", "test"))
    nls_tr <- fmt(get_mean(TRUE, "oracle-nl",     "train")); nls_te <- fmt(get_mean(TRUE, "oracle-nl",     "test"))
    ls_tr  <- fmt(get_mean(TRUE, "lasso",         "train")); ls_te  <- fmt(get_mean(TRUE, "lasso",         "test"))

    out <- paste0(out, "\n", glue::glue(
"================= Stratified ==================
(Oracle linear)     Train: {os_tr} |  Test: {os_te}
(Oracle nonlinear)  Train: {nls_tr} |  Test: {nls_te}
(Lasso)             Train: {ls_tr} |  Test: {ls_te}"
    ))
  }
  paste0(glue::glue("Results for N={n_train}:\n\n"), out, "\n\n")
}


#  ---- Try print out ------------------------

# * ---- linear moderate ----
result_files <- list.files(here::here('results', 'models', 'BE10Asso00_normTC', 'linear-moderate','cox_all'), full.names = T)
results <- read.csv(result_files[[1]])

cat(
  unlist(c(100, 200, 500) |> 
    purrr::map(summarize_block, results = results) 
))

# * ---- nl quadratic ----
result_files <- list.files(here::here('results', 'models', 'BE10Asso00_normTC', 'nl-quadratic','cox_all'), full.names = T)
results <- read.csv(result_files[[1]])

cat(unlist(
  c(100, 200, 500) |> 
    purrr::map(summarize_block, results = results) 
))

