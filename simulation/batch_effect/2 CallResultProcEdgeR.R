
source("~/dl-survival-miRNA/simulation/batch_effect/resultproc.R")

################## large effect ###################

### no train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he0000bs20norm00tsort00", output_name="edgeRV2he0000bs20norm00tsort00")
# 
# ### no train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=1, exclude=100, addon="he0000bs20norm11tsort00", output_name="edgeRV2he0000bs20norm11tsort00")
# 
# ## no train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he0000bs20norm22tsort00", output_name="edgeRV2he0000bs20norm22tsort00")
# 
# ### no train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he0000bs20norm33tsort00", output_name="edgeRV2he0000bs20norm33tsort00")
# 
# ### no train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he0000bs20norm44tsort00", output_name="edgeRV2he0000bs20norm44tsort00")



### train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm00tsort00", output_name="edgeRV2he4400bs20norm00tsort00")
# 
# ### train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm11tsort00", output_name="edgeRV2he4400bs20norm11tsort00")
# 
# ## train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm22tsort00", output_name="edgeRV2he4400bs20norm22tsort00")
# 
# ### train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm33tsort00", output_name="edgeRV2he4400bs20norm33tsort00")
# 
# ### train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm44tsort00", output_name="edgeRV2he4400bs20norm44tsort00")



### train HE, no test HE, sort train, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm00tsort0.050", output_name="edgeRV2he4400bs20norm00tsort0.050")
# 
# ### train HE, no test HE, sort train, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm11tsort0.050", output_name="edgeRV2he4400bs20norm11tsort0.050")
# 
# ## train HE, no test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm22tsort0.050", output_name="edgeRV2he4400bs20norm22tsort0.050")
# 
# ### train HE, no test HE, sort train, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm33tsort0.050", output_name="edgeRV2he4400bs20norm33tsort0.050")
# 
# ### train HE, no test HE, sort train, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4400bs20norm44tsort0.050", output_name="edgeRV2he4400bs20norm44tsort0.050")



### train HE, test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=80, exclude=100, addon="he4444bs20norm00tsort00", output_name="edgeRV2he4444bs20norm00")
# 
# ### train HE, test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4444bs20norm11tsort00", output_name="edgeRV2he4444bs20norm11")
# 
# ## train HE, test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4444bs20norm22tsort00", output_name="edgeRV2he4444bs20norm22")
# 
# ### train HE, test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4444bs20norm33tsort00", output_name="edgeRV2he4444bs20norm33")
# 
# ### train HE, test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=20, exclude=100, addon="he4444bs20norm44tsort00", output_name="edgeRV2he4444bs20norm44")



### train HE, test HE, sort train, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm00tsort0.050", output_name="edgeRV2he4444bs20norm00tsort0.050")
# 
# ### train HE, test HE, sort train, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm11tsort0.050", output_name="edgeRV2he4444bs20norm11tsort0.050")

# ## train HE, test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm22tsort0.050", output_name="edgeRV2he4444bs20norm22tsort0.050")
# 
# ### train HE, test HE, sort train, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm33tsort0.050", output_name="edgeRV2he4444bs20norm33tsort0.050")
# 
# ### train HE, test HE, sort train, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm44tsort0.050", output_name="edgeRV2he4444bs20norm44tsort0.050")



### train HE, test HE, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm00tsort00.05", output_name="edgeRV2he4444bs20norm00tsort00.05")
# 
# ### train HE, test HE, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm11tsort00.05", output_name="edgeRV2he4444bs20norm11tsort00.05")

# ## train HE, test HE, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm22tsort00.05", output_name="edgeRV2he4444bs20norm22tsort00.05")
# 
# ### train HE, test HE, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm33tsort00.05", output_name="edgeRV2he4444bs20norm33tsort00.05")
# 
# ### train HE, test HE, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm44tsort00.05", output_name="edgeRV2he4444bs20norm44tsort00.05")



### train HE, test HE, sort train, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm00tsort0.050.05", output_name="edgeRV2he4444bs20norm00tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm11tsort0.050.05", output_name="edgeRV2he4444bs20norm11tsort0.050.05")

## train HE, test HE, sort train, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm22tsort0.050.05", output_name="edgeRV2he4444bs20norm22tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm33tsort0.050.05", output_name="edgeRV2he4444bs20norm33tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm44tsort0.050.05", output_name="edgeRV2he4444bs20norm44tsort0.050.05")



### train HE, test HE, sort train, inverse sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm00tsort0.05-0.05", output_name="edgeRV2he4444bs20norm00tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm11tsort0.05-0.05", output_name="edgeRV2he4444bs20norm11tsort0.05-0.05")

# ## train HE, test HE, sort train, inverse sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm22tsort0.05-0.05", output_name="edgeRV2he4444bs20norm22tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm33tsort0.05-0.05", output_name="edgeRV2he4444bs20norm33tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he4444bs20norm44tsort0.05-0.05", output_name="edgeRV2he4444bs20norm44tsort0.05-0.05")




















###################################### Large gene effects, mean-level batch effects only ###################################

# ### no train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm00tsort00", output_name="edgeRV2he0000bs20norm00tsort00")
# 
# ### no train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm11tsort00", output_name="edgeRV2he0000bs20norm11tsort00")
# 
# ## no train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm22tsort00", output_name="edgeRV2he0000bs20norm22tsort00")
# 
# ### no train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm33tsort00", output_name="edgeRV2he0000bs20norm33tsort00")
# 
# ### no train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm44tsort00", output_name="edgeRV2he0000bs20norm44tsort00")
# 
# 
# 
# ### train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm00tsort00", output_name="edgeRV2he2000bs20norm00tsort00")
# 
# ### train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm11tsort00", output_name="edgeRV2he2000bs20norm11tsort00")
# 
# ## train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm22tsort00", output_name="edgeRV2he2000bs20norm22tsort00")
# 
# ### train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm33tsort00", output_name="edgeRV2he2000bs20norm33tsort00")
# 
# ### train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm44tsort00", output_name="edgeRV2he2000bs20norm44tsort00")



# ### train HE, no test HE, sort train, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm00tsort0.050", output_name="edgeRV2he2000bs20norm00tsort0.050")
# 
# ### train HE, no test HE, sort train, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm11tsort0.050", output_name="edgeRV2he2000bs20norm11tsort0.050")
# 
# ## train HE, no test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm22tsort0.050", output_name="edgeRV2he2000bs20norm22tsort0.050")
# 
# ### train HE, no test HE, sort train, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm33tsort0.050", output_name="edgeRV2he2000bs20norm33tsort0.050")
# 
# ### train HE, no test HE, sort train, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2000bs20norm44tsort0.050", output_name="edgeRV2he2000bs20norm44tsort0.050")
# 
# 
# 
# ### train HE, test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2020bs20norm00tsort00", output_name="edgeRV2he2020bs20norm00tsort00")
# 
# ### train HE, test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2020bs20norm11tsort00", output_name="edgeRV2he2020bs20norm11tsort00")
# 
# ## train HE, test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2020bs20norm22tsort00", output_name="edgeRV2he2020bs20norm22tsort00")
# 
# ### train HE, test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2020bs20norm33tsort00", output_name="edgeRV2he2020bs20norm33tsort00")
# 
# ### train HE, test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2020bs20norm44tsort00", output_name="edgeRV2he2020bs20norm44tsort00")

# 
# 
# ### train HE, test HE, sort train, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm00tsort0.050", output_name="edgeRV2he2020bs20norm00tsort0.050")
# 
# ### train HE, test HE, sort train, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm11tsort0.050", output_name="edgeRV2he2020bs20norm11tsort0.050")
# 
# ## train HE, test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm22tsort0.050", output_name="edgeRV2he2020bs20norm22tsort0.050")
# 
# ### train HE, test HE, sort train, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm33tsort0.050", output_name="edgeRV2he2020bs20norm33tsort0.050")
# 
# ### train HE, test HE, sort train, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm44tsort0.050", output_name="edgeRV2he2020bs20norm44tsort0.050")
# 
# 
# 
# ### train HE, test HE, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm00tsort00.05", output_name="edgeRV2he2020bs20norm00tsort00.05")
# 
# ### train HE, test HE, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm11tsort00.05", output_name="edgeRV2he2020bs20norm11tsort00.05")
# 
# ## train HE, test HE, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm22tsort00.05", output_name="edgeRV2he2020bs20norm22tsort00.05")
# 
# ### train HE, test HE, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm33tsort00.05", output_name="edgeRV2he2020bs20norm33tsort00.05")
# 
# ### train HE, test HE, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm44tsort00.05", output_name="edgeRV2he2020bs20norm44tsort00.05")
# 
# 
# 
# ### train HE, test HE, sort train, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm00tsort0.050.05", output_name="edgeRV2he2020bs20norm00tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm11tsort0.050.05", output_name="edgeRV2he2020bs20norm11tsort0.050.05")
# 
# ## train HE, test HE, sort train, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm22tsort0.050.05", output_name="edgeRV2he2020bs20norm22tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm33tsort0.050.05", output_name="edgeRV2he2020bs20norm33tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm44tsort0.050.05", output_name="edgeRV2he2020bs20norm44tsort0.050.05")
# 
# 
# 
# ### train HE, test HE, sort train, inverse sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm00tsort0.05-0.05", output_name="edgeRV2he2020bs20norm00tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm11tsort0.05-0.05", output_name="edgeRV2he2020bs20norm11tsort0.05-0.05")
# 
# ## train HE, test HE, sort train, inverse sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm22tsort0.05-0.05", output_name="edgeRV2he2020bs20norm22tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm33tsort0.05-0.05", output_name="edgeRV2he2020bs20norm33tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2020bs20norm44tsort0.05-0.05", output_name="edgeRV2he2020bs20norm44tsort0.05-0.05")










######################################### Large gene effects, observation-level batch effects only ######################################

# ### no train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm00tsort00", output_name="edgeRV2he0000bs20norm00tsort00")
# 
# ### no train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm11tsort00", output_name="edgeRV2he0000bs20norm11tsort00")
# 
# ## no train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm22tsort00", output_name="edgeRV2he0000bs20norm22tsort00")
# 
# ### no train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm33tsort00", output_name="edgeRV2he0000bs20norm33tsort00")
# 
# ### no train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm44tsort00", output_name="edgeRV2he0000bs20norm44tsort00")
# 
# 
# 
# ### train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm00tsort00", output_name="edgeRV2he0200bs20norm00tsort00")
# 
# ### train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm11tsort00", output_name="edgeRV2he0200bs20norm11tsort00")
# 
# ## train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm22tsort00", output_name="edgeRV2he0200bs20norm22tsort00")
# 
# ### train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm33tsort00", output_name="edgeRV2he0200bs20norm33tsort00")
# 
# ### train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm44tsort00", output_name="edgeRV2he0200bs20norm44tsort00")



# ### train HE, no test HE, sort train, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm00tsort0.050", output_name="edgeRV2he0200bs20norm00tsort0.050")
# 
# ### train HE, no test HE, sort train, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm11tsort0.050", output_name="edgeRV2he0200bs20norm11tsort0.050")
# 
# ## train HE, no test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm22tsort0.050", output_name="edgeRV2he0200bs20norm22tsort0.050")
# 
# ### train HE, no test HE, sort train, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm33tsort0.050", output_name="edgeRV2he0200bs20norm33tsort0.050")
# 
# ### train HE, no test HE, sort train, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0200bs20norm44tsort0.050", output_name="edgeRV2he0200bs20norm44tsort0.050")
# 
# 
# 
# ### train HE, test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0202bs20norm00tsort00", output_name="edgeRV2he0202bs20norm00tsort00")
# 
# ### train HE, test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0202bs20norm11tsort00", output_name="edgeRV2he0202bs20norm11tsort00")
# 
# ## train HE, test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0202bs20norm22tsort00", output_name="edgeRV2he0202bs20norm22tsort00")
# 
# ### train HE, test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0202bs20norm33tsort00", output_name="edgeRV2he0202bs20norm33tsort00")
# 
# ### train HE, test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0202bs20norm44tsort00", output_name="edgeRV2he0202bs20norm44tsort00")


# 
# ### train HE, test HE, sort train, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm00tsort0.050", output_name="edgeRV2he0202bs20norm00tsort0.050")
# 
# ### train HE, test HE, sort train, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm11tsort0.050", output_name="edgeRV2he0202bs20norm11tsort0.050")
# 
# ## train HE, test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm22tsort0.050", output_name="edgeRV2he0202bs20norm22tsort0.050")
# 
# ### train HE, test HE, sort train, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm33tsort0.050", output_name="edgeRV2he0202bs20norm33tsort0.050")
# 
# ### train HE, test HE, sort train, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm44tsort0.050", output_name="edgeRV2he0202bs20norm44tsort0.050")
# 
# 
# 
# ### train HE, test HE, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm00tsort00.05", output_name="edgeRV2he0202bs20norm00tsort00.05")
# 
# ### train HE, test HE, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm11tsort00.05", output_name="edgeRV2he0202bs20norm11tsort00.05")
# 
# ## train HE, test HE, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm22tsort00.05", output_name="edgeRV2he0202bs20norm22tsort00.05")
# 
# ### train HE, test HE, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm33tsort00.05", output_name="edgeRV2he0202bs20norm33tsort00.05")
# 
# ### train HE, test HE, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm44tsort00.05", output_name="edgeRV2he0202bs20norm44tsort00.05")
# 
# 
# 
# ### train HE, test HE, sort train, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm00tsort0.050.05", output_name="edgeRV2he0202bs20norm00tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm11tsort0.050.05", output_name="edgeRV2he0202bs20norm11tsort0.050.05")
# 
# ## train HE, test HE, sort train, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm22tsort0.050.05", output_name="edgeRV2he0202bs20norm22tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm33tsort0.050.05", output_name="edgeRV2he0202bs20norm33tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm44tsort0.050.05", output_name="edgeRV2he0202bs20norm44tsort0.050.05")
# 
# 
# 
# ### train HE, test HE, sort train, inverse sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm00tsort0.05-0.05", output_name="edgeRV2he0202bs20norm00tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm11tsort0.05-0.05", output_name="edgeRV2he0202bs20norm11tsort0.05-0.05")
# 
# ## train HE, test HE, sort train, inverse sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm22tsort0.05-0.05", output_name="edgeRV2he0202bs20norm22tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm33tsort0.05-0.05", output_name="edgeRV2he0202bs20norm33tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he0202bs20norm44tsort0.05-0.05", output_name="edgeRV2he0202bs20norm44tsort0.05-0.05")









######################################### Large gene effects, both mean-level and observation-level batch effects ######################################

# ### no train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm00tsort00", output_name="edgeRV2he0000bs20norm00tsort00")
# 
# ### no train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm11tsort00", output_name="edgeRV2he0000bs20norm11tsort00")
# 
# ## no train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm22tsort00", output_name="edgeRV2he0000bs20norm22tsort00")
# 
# ### no train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm33tsort00", output_name="edgeRV2he0000bs20norm33tsort00")
# 
# ### no train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he0000bs20norm44tsort00", output_name="edgeRV2he0000bs20norm44tsort00")
# 
# 
# 
# ### train HE, no test HE, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm00tsort00", output_name="edgeRV2he2200bs20norm00tsort00")
# 
# ### train HE, no test HE, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm11tsort00", output_name="edgeRV2he2200bs20norm11tsort00")
# 
# ## train HE, no test HE, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm22tsort00", output_name="edgeRV2he2200bs20norm22tsort00")
# 
# ### train HE, no test HE, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm33tsort00", output_name="edgeRV2he2200bs20norm33tsort00")
# 
# ### train HE, no test HE, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm44tsort00", output_name="edgeRV2he2200bs20norm44tsort00")

### train HE, no test HE, sort train, no norm, batch size 20, even larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm00tsort0.050", output_name="edgeRV2he2200bs20norm00tsort0.050")

### train HE, no test HE, sort train, TC norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm11tsort0.050", output_name="edgeRV2he2200bs20norm11tsort0.050")

## train HE, no test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm22tsort0.050", output_name="edgeRV2he2200bs20norm22tsort0.050")

### train HE, no test HE, sort train, TMM norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm33tsort0.050", output_name="edgeRV2he2200bs20norm33tsort0.050")

### train HE, no test HE, sort train, quantile norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2200bs20norm44tsort0.050", output_name="edgeRV2he2200bs20norm44tsort0.050")



### train HE, test HE, no norm, batch size 20, even larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2222bs20norm00tsort00", output_name="edgeRV2he2222bs20norm00tsort00")

### train HE, test HE, TC norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2222bs20norm11tsort00", output_name="edgeRV2he2222bs20norm11tsort00")

## train HE, test HE, upper quantile norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2222bs20norm22tsort00", output_name="edgeRV2he2222bs20norm22tsort00")

### train HE, test HE, TMM norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2222bs20norm33tsort00", output_name="edgeRV2he2222bs20norm33tsort00")

### train HE, test HE, quantile norm, batch size 20, larger batch effect ###
resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
           stratify=1, cmbt=0, nbatch=20, exclude=100, addon="he2222bs20norm44tsort00", output_name="edgeRV2he2222bs20norm44tsort00")



# ### train HE, test HE, sort train, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm00tsort0.050", output_name="edgeRV2he2222bs20norm00tsort0.050")
# 
# ### train HE, test HE, sort train, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm11tsort0.050", output_name="edgeRV2he2222bs20norm11tsort0.050")
# 
# ## train HE, test HE, sort train, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm22tsort0.050", output_name="edgeRV2he2222bs20norm22tsort0.050")
# 
# ### train HE, test HE, sort train, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm33tsort0.050", output_name="edgeRV2he2222bs20norm33tsort0.050")
# 
# ### train HE, test HE, sort train, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm44tsort0.050", output_name="edgeRV2he2222bs20norm44tsort0.050")
# 
# 
# 
# ### train HE, test HE, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm00tsort00.05", output_name="edgeRV2he2222bs20norm00tsort00.05")
# 
# ### train HE, test HE, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm11tsort00.05", output_name="edgeRV2he2222bs20norm11tsort00.05")
# 
# ## train HE, test HE, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm22tsort00.05", output_name="edgeRV2he2222bs20norm22tsort00.05")
# 
# ### train HE, test HE, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm33tsort00.05", output_name="edgeRV2he2222bs20norm33tsort00.05")
# 
# ### train HE, test HE, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm44tsort00.05", output_name="edgeRV2he2222bs20norm44tsort00.05")
# 
# 
# 
# ### train HE, test HE, sort train, sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm00tsort0.050.05", output_name="edgeRV2he2222bs20norm00tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm11tsort0.050.05", output_name="edgeRV2he2222bs20norm11tsort0.050.05")
# 
# ## train HE, test HE, sort train, sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm22tsort0.050.05", output_name="edgeRV2he2222bs20norm22tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm33tsort0.050.05", output_name="edgeRV2he2222bs20norm33tsort0.050.05")
# 
# ### train HE, test HE, sort train, sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm44tsort0.050.05", output_name="edgeRV2he2222bs20norm44tsort0.050.05")
# 
# 
# 
# ### train HE, test HE, sort train, inverse sort test, no norm, batch size 20, even larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm00tsort0.05-0.05", output_name="edgeRV2he2222bs20norm00tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TC norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm11tsort0.05-0.05", output_name="edgeRV2he2222bs20norm11tsort0.05-0.05")
# 
# ## train HE, test HE, sort train, inverse sort test, upper quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm22tsort0.05-0.05", output_name="edgeRV2he2222bs20norm22tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, TMM norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm33tsort0.05-0.05", output_name="edgeRV2he2222bs20norm33tsort0.05-0.05")
# 
# ### train HE, test HE, sort train, inverse sort test, quantile norm, batch size 20, larger batch effect ###
# resultproc(nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444), b0value=c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9), p=462,
#            stratify=1, cmbt=1, nbatch=40, exclude=100, addon="he2222bs20norm44tsort0.05-0.05", output_name="edgeRV2he2222bs20norm44tsort0.05-0.05")

