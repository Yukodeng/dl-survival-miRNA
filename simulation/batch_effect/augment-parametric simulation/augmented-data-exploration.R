
setwd("~/dl-survival-miRNA/simulation/batch_effect/with_without_batch_maf/")
data = mxf
data_batch = read.csv('batch_10_epochES_batch01_maf_generated.csv')

hist(data$hsa.let.7a..2..1)
hist(data_batch$hsa.let.7a..2..1)
