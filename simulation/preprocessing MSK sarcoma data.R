
############ use edgeR to estimate mu and phi for each gene using original MSK dataset (n=27+27)
library(edgeR)

# setwd("Augmented sequencing data analysis")
setwd("~/dl-survival-miRNA/simulation")

# mirna.names = colnames(read.csv("MSKpair_MXF_merge.csv")[2:2067])
# mirna.names = as.data.frame(mirna.names)
# write.csv(mirna.names, "../data/MSKpair_miRNA_names.csv")

mxf.rwd=read.csv("MSKpair_MXF_merge.csv")[,2:1034]
pmfh.rwd=read.csv("MSKpair_PMFH_merge.csv")[,2:1034]
# working.rwd=rbind(mxf.rwd, pmfh.rwd)
workingt.rwd=t(rbind(mxf.rwd, pmfh.rwd))
cpm.rwd=cpm(workingt.rwd)

keep=(rowSums(cpm.rwd>2)>=5) ## results in 667 genes
workingt.rwd=workingt.rwd[keep,]
n=ncol(workingt.rwd)
p=nrow(workingt.rwd)
cpm.rwd=cpm.rwd[keep,]
cor.mat=cor(t(cpm))
design=rep(1,ncol(workingt))
# 
# disps=estimateDisp(workingt, design)
# disps.gene=disps$tagwise.dispersion
# fit=glmFit(workingt, dispersion=disps.gene)  ## by default, design matrix is a vector of ones, no offset, library size is log(colSums(workingt))
# betas.gene=as.numeric(fit$coefficients)
# libsizes=colSums(workingt)
# 
# write.csv(betas.gene, "betas.gene.csv", row.names=F)
# write.csv(libsizes, "libsizes.csv", row.names=F)
# write.csv(disps.gene, "disps.gene.csv", row.names=F)
# write.csv(cor.mat, "cor.mat.csv", row.names=F)

## select 10 true genes associated with survival time
tc=rowSums(workingt.rwd)
mid.tc=(tc>median(tc) & tc<quantile(tc, 0.75)) # select genes with total count in the third quartile
cv=apply(workingt.rwd, 1, sd)/rowMeans(workingt.rwd)
mid.cv=(cv>median(cv)) # select genes with coefficient of variation (SD/mean) in the upper half
geneid=row.names(workingt.rwd)
candidate.geneid=geneid[mid.tc & mid.cv]
set.seed(1234)
selected.geneid=sample(candidate.geneid, 10, replace=F)
nonzero_position=which(geneid %in% selected.geneid)  ## selected gene position c(71, 72, 101, 123, 136, 166, 215, 354, 577, 640)

selected.geneid
nonzero_position

# selected.gene=workingt[nonzero_position,]
# rowMeans(selected.gene)
# sel.sd=apply(selected.gene, 1, sd)
# sel.cor=cor(t(selected.gene))
# summary(sel.cor[sel.cor!=1])

selected.gene.cpm=cpm.rwd[nonzero_position,]
sel.sd.cpm=apply(log2(selected.gene.cpm+0.5), 1, sd)


beta0=rep(c(1,-1),5)*round(5/sel.sd.cpm, 1)  ## c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9)
# write.csv(data.frame(beta0), paste0("../data/beta0.csv"))

beta0


selected.gene.cpm.log = log(selected.gene.cpm+0.5)
View(scale(selected.gene.cpm.log))
x.mean = rowMeans(scale(selected.gene.cpm.log))
write.csv(data.frame(x.mean), 'selected.gene.mean.scaled.csv')
