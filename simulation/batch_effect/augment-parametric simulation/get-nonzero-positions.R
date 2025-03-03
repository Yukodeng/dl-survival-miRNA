setwd('~/dl-survival-miRNA/simulation/batch_effect/')

# load selected genes
geneid=names(read.csv("cor.mat.462.csv"))  ## estimated correlation matrix among genes
nonzero_position=c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444)
sel_genes = geneid[nonzero_position] 

# new set of genes
new.geneid = read.csv(file.path('augment-parametric simulation', 'gene.names.538.csv'))$gene

# retrive the new positions of the target genes
new.nonzero_position = which(new.geneid %in% sel_genes)

all(new.geneid[new.nonzero_position] == sel_genes) # make sure this is true

print(new.nonzero_position)
