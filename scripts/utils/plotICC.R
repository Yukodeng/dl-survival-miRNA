library(lme4)
library(ggplot2)
library(dplyr)

compute_icc <- function(y, batch.id) {
  if (length(table(y)) == 1) return(0)
  fit <- suppressMessages(lmer(y ~ 1 + (1 | batch.id)))
  vc <- VarCorr(fit)
  sigma_b <- as.numeric(vc$batch.id)       # between-group variance
  sigma_w <- attr(vc, "sc")^2              # residual variance
  return(sigma_b / (sigma_b + sigma_w))
}

# X: n x p matrix
# batch.id: length n
compute_icc_matrix = function(X, batch.id){
  batch.id <- factor(batch.id)
  return(apply(X, 2, compute_icc, batch.id = batch.id))
}

## MARK
# real = read.csv("./data/MSK sarcoma with and without 15 batch counts 01092026_CVAE_0.01_1000perbatch.csv")

## 02052026 read in data
real.clean.log = read.csv("par_clean.csv")
real.batch.log = read.csv("par_batch.csv")
aug.clean.log = read.csv("aug_clean.csv")
aug.batch.log = read.csv("aug_batch.csv")

batch.num = 15             # number of batches
sample.per.batch = 1000    # number of samples per batch


# clean.cols = filter.idx
# batch.cols = filter.idx + 1033
# real.clean = real[, clean.cols]
# real.batch = real[, batch.cols]

# real.batch.log = (real.batch + 1) %>%
#   log2() %>%
#   as.matrix()
# real.clean.log = (real.clean + 1) %>%
#   log2() %>%
#   as.matrix()

## parametric  ---------------------------------

# pure.batch = real.batch.log - real.clean.log
# sd.clean = apply(real.clean.log, 2, sd)
# sd.pure.batch = apply(pure.batch, 2, sd)
# mean.clean = apply(real.clean.log, 2, mean)
# mean.batch = apply(real.batch.log, 2, mean)
batch.id = rep(1:batch.num, each = 54)

icc.batch = compute_icc_matrix(real.batch, batch.id)
# na.idx = is.na(icc.batch)
mean.batch = colMeans(real.batch.log)
sd.batch = sapply(1:ncol(real.batch.log), function(x) sd(real.batch.log[, x]))
df.icc <- data.frame(
  mean_expr = mean.batch,
  icc = icc.batch
)
p = ggplot(df.icc, aes(x = mean_expr, y = icc)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    x = "Log2 mean expression",
    y = "ICC (batch)",
    title = "Feature-wise ICC vs Mean Expression (Parametric)"
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    color = "blue",
    linewidth = 1
  ) +
  theme_classic()
print(p)
# ggsave("./plot/ICC-Mean.png", plot = p, width = 6, height = 5, dpi = 300)



## augmented -------------------------------------

pure.batch = aug.batch.log - aug.clean.log
sd.clean = apply(aug.clean.log, 2, sd)
sd.pure.batch = apply(pure.batch, 2, sd)
mean.clean = apply(aug.clean.log, 2, mean)
mean.batch = apply(aug.batch.log, 2, mean)

batch.id = rep(1:batch.num, each = sample.per.batch)

icc.batch = compute_icc_matrix(aug.batch, batch.id)
# na.idx = is.na(icc.batch)
mean.batch = colMeans(aug.batch.log)
sd.batch = sapply(1:ncol(aug.batch.log), function(x) sd(aug.batch.log[, x]))
df.icc <- data.frame(
  mean_expr = mean.batch,
  icc = icc.batch
)
p = ggplot(df.icc, aes(x = mean_expr, y = icc)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    x = "Log2 mean expression",
    y = "ICC (batch)",
    title = "Feature-wise ICC vs Mean Expression (Augmented)"
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    color = "blue",
    linewidth = 1
  ) +
  theme_classic()
print(p)
ggsave("./plot/ICC-Mean.png", plot = p, width = 6, height = 5, dpi = 300)