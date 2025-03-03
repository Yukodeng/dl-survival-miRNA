
source("~/dl-survival-miRNA/simulation/batch_effect/resultproc.R")

###################################### Large gene effects, mean-level batch effects only ###################################

### no train HE, no test HE, no norm, batch size 20, even larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm00tsort00", output_name="edgeRV2he0000bs20norm00tsort00")

### no train HE, no test HE, TC norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm11tsort00", output_name="edgeRV2he0000bs20norm11tsort00")

## no train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm22tsort00", output_name="edgeRV2he0000bs20norm22tsort00")
