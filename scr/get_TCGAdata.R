# BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

project.name = 'TCGA-SKCM'
# getProjectSummary("TCGA-SKCM")

# get_TCGAdata <- function(project.name, directory) {
#   
# }

# Query for miRNA expression data ----
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

# get miRNA expression matrix
skcm.miRNA.counts <-miRNA.data %>%
  select(colnames(miRNA.data)[grepl("reads_per_million|miRNA_ID", colnames(miRNA.data))]) %>% 
  tibble::column_to_rownames("miRNA_ID") 
colnames(skcm.miRNA.counts) <- substr(gsub(".*TCGA", "TCGA", colnames(skcm.miRNA.counts)), 0, 12)


# Get TCGA project clinical data ----
clinical <- GDCquery_clinic(
  project = project.name,
  type = "clinical"
)


# Save miRNA expression data to a CSV file
write.csv(skcm.miRNA.counts, 
          file.path("data", "batch", paste0(project.name, "_miRNA_counts.csv"))
)
# Save clinical survival data to a CSV file
write.csv(clinical, 
          file.path("data", "batch", paste0(project.name, "_clinical_data.csv")), 
          row.names=FALSE)
