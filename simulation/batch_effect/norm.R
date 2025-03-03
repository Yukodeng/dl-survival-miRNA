norm.TC <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  dat.DGE <- edgeR::DGEList(
    counts = matrix(raw, ncol = length(groups)),
    group = factor(groups),
    genes = rownames(raw)
  )
  scalingFactor <- dat.DGE$samples$lib.size / 1e6
  dataNormalized <- t(t(raw) / scalingFactor)
  # scalingFactor <- dat.DGE$samples$lib.size / mean(dat.DGE$samples$lib.size)
  # dataNormalized <- round(t(t(raw)/scalingFactor))
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor
  ))
}

norm.UQ <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  dat.DGE <- edgeR::DGEList(
    counts = matrix(raw, ncol = length(groups)),
    group = factor(groups),
    genes = rownames(raw)
  )
  q.factor <- apply(dat.DGE$counts, 2, function(x) stats::quantile(x[x != 0], probs = 0.75))
  # scalingFactor <- q.factor / 1e6
  scalingFactor <- q.factor / 1e3
  dataNormalized <- t(t(raw) / scalingFactor)
  # scalingFactor <- q.factor/mean(dat.DGE$samples$lib.size)
  # dataNormalized <- round(t(t(raw)/scalingFactor))
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor
  ))
}

norm.TMM <- function(raw, groups = rep(1, ncol(raw))) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package \"edgeR\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  dat.DGE <- edgeR::DGEList(
    counts = matrix(raw, ncol = length(groups)),
    group = factor(groups),
    genes = rownames(raw)
  )
  d <- edgeR::calcNormFactors(dat.DGE, method = "TMM")
  scalingFactor <- d$samples$norm.factors * d$samples$lib.size / 1e6
  dataNormalized <- t(t(raw) / scalingFactor)
  # scalingFactor <- d$samples$norm.factors * d$samples$lib.size /
  #   mean(dat.DGE$samples$lib.size)
  # dataNormalized <- round(t(t(raw)/scalingFactor))
  return(list(
    dataNormalized = dataNormalized,
    scalingFactor = scalingFactor
  ))
}

norm.PoissonSeq <- function(raw) {
  invisible(capture.output(scaling.factor <- PS.Est.Depth(raw)))
  dat.normed <- t(t(raw)/scaling.factor)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}
