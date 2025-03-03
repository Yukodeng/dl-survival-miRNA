source("~/dl-survival-miRNA/simulation/batch_effect/resultproc.R")
p=538

###################################### Large gene effects, mean-level batch effects only ###################################
### no train HE, no test HE, no norm, batch size 1, even larger batch effect ###
######## NOTE: different folder/file names is the only exception! Change to align with the rest later on
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he00bs20norm00tsort00-aug", output_name="edgeRV2he00bs20norm00tsort00-aug")

### no train HE, no test HE, TC norm, batch size 1, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he0000bs20norm11tsort00-aug", output_name="edgeRV2he0000bs20norm11tsort00-aug")

### no train HE, no test HE, upper quantile norm, batch size 1, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he0000bs20norm22tsort00-aug", output_name="edgeRV2he0000bs20norm22tsort00-aug")

### no train HE, no test HE, TMM norm, batch size 1, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he0000bs20norm33tsort00-aug", output_name="edgeRV2he0000bs20norm33tsort00-aug")

### no train HE, no test HE, quantile norm, batch size 1, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he0000bs20norm44tsort00-aug", output_name="edgeRV2he0000bs20norm44tsort00-aug")



#########################################################################################
### train HE, no test HE, sort train, no norm, batch size 20, even larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he2000bs20norm00tsort00-aug", output_name="edgeRV2he2000bs20norm00tsort00-aug")

### train HE, no test HE, sort train, TC norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he2000bs20norm11tsort00-aug", output_name="edgeRV2he2000bs20norm11tsort00-aug")

## train HE, no test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he2000bs20norm22tsort00-aug", output_name="edgeRV2he2000bs20norm22tsort00-aug")

### train HE, no test HE, sort train, TMM norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he2000bs20norm33tsort0-aug", output_name="edgeRV2he2000bs20norm33tsort00-aug")

### train HE, no test HE, sort train, quantile norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=p,
           stratify=1, cmbt=0, nbatch=1, exclude=100, addon="he2000bs20norm44tsort00-aug", output_name="edgeRV2he2000bs20norm44tsort00-aug")
