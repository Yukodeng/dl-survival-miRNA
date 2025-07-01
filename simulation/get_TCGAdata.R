# BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(preprocessCore)
library(dplyr)


get_TCGAdata <- function(project.name, save.data = T, ignore.exist=F) {
  
  if ( length(list.files(file.path("data"), paste0(project.name,"_miRNA_counts"))) != 0 &
       length(list.files(file.path("data"), paste0(project.name,"_clinical"))) != 0 &
       !ignore.exist
       ) {
    print("Files already exist in the directory.")
    return (NULL)
  }
  
  # Get miRNA expression data ------------------
  miRNA.query <- GDCquery(
    project = project.name,                  # TCGA project
    data.category = "Transcriptome Profiling",  # Data category
    data.type = "miRNA Expression Quantification", # Data type
    workflow.type = "BCGSC miRNA Profiling", # Workflow type
    experimental.strategy = "miRNA-Seq"      # Experimental strategy
  )
  # download from GDC and prepare data
  GDCdownload(query = miRNA.query)
  miRNA.data <- GDCprepare(query = miRNA.query)
  
 
  ## Preprocess miRNA expression count data ----
  # 1. data cleaning
  miRNA.counts <-miRNA.data %>%
    select(colnames(miRNA.data)[grepl("reads_per_million|miRNA_ID", colnames(miRNA.data))]) %>% 
    tibble::column_to_rownames("miRNA_ID") 
  colnames(miRNA.counts) <- substr(gsub(".*TCGA", "TCGA", colnames(miRNA.counts)), 0, 12)
  
  # 2. quantile normalization
  miRNA.norm = normalize.quantiles(apply(miRNA.counts, 2, as.numeric), keep.names = T )
  rownames(miRNA.norm) = rownames(miRNA.counts)
  
  # 3. log2-transformation
  miRNA.norm = log2(miRNA.norm +1)
  
  
  # Get clinical data ----------------------------
  clinical <- GDCquery_clinic(
    project = project.name,
    type = "clinical"
  )

  
  # Save Data ------------------------------------
  if (save.data) {
    ## Save miRNA expression data to a CSV file
    write.csv(miRNA.norm, 
              file.path("data", "batch", paste0(project.name, "_miRNA_counts.csv"))
    )
 
    ## Save clinical survival data to a CSV file
    write.csv(clinical, 
              file.path("data", "batch", paste0(project.name, "_clinical_data.csv")), 
              row.names=FALSE
    )
  }
  
  out = list(miRNA.norm, clinical)
  return(out)
}


# Example: Download TCGA-SKCM data ----------
# project.name = 'TCGA-SKCM'
# # getProjectSummary("TCGA-SKCM")
# get_TCGAdata(project.name = project.name, save.data = F)


## LAML data
get_TCGAdata(project.name = 'TCGA-LAML', save.data = T, ignore.exist = T)


### do it manually ======= 
# project.name = 'TCGA-LAML'
# # Get miRNA expression data ------------------
# miRNA.query <- GDCquery(
#   project = project.name,                  # TCGA project
#   data.category = "Transcriptome Profiling",  # Data category
#   data.type = "miRNA Expression Quantification", # Data type
#   workflow.type = "BCGSC miRNA Profiling", # Workflow type
#   experimental.strategy = "miRNA-Seq"      # Experimental strategy
# )
# # download from GDC and prepare data
# GDCdownload(query = miRNA.query)
# miRNA.data <- GDCprepare(query = miRNA.query)
# 
# 
# ## Preprocess miRNA expression count data ----
# # 1. data cleaning
# miRNA.counts <-miRNA.data %>%
#   select(colnames(miRNA.data)[grepl("reads_per_million|miRNA_ID", colnames(miRNA.data))]) %>% 
#   tibble::column_to_rownames("miRNA_ID") 
# colnames(miRNA.counts) <- substr(gsub(".*TCGA", "TCGA", colnames(miRNA.counts)), 0, 12)
# 
# # 2. quantile normalization
# miRNA.norm = normalize.quantiles(apply(miRNA.counts, 2, as.numeric), keep.names = T )
# rownames(miRNA.norm) = rownames(miRNA.counts)
# 
# # 3. log2-transformation
# miRNA.norm = log2(miRNA.norm +1)
# 
# 
# # Get clinical data ----------------------------
# clinical <- GDCquery_clinic(
#   project = project.name,
#   type = "clinical"
# )
# 
# 
# # Save Data ------------------------------------
# ## Save miRNA expression data to a CSV file
# write.csv(miRNA.counts, 
#           file.path("data", "batch", paste0(project.name, "_miRNA_counts.csv")),
# )
#   
# ## Save clinical survival data to a CSV file
# write.csv(clinical, 
#           file.path("data", "batch", paste0(project.name, "_clinical_data.csv")), 
#           row.names=FALSE
# )