# Load packages
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(SummarizedExperiment)


# Set directory
setwd("D:/R_projects/mutual_exc/tcga_data/")
source(file = "scripts/maf2matrix.R")

# TCGA
TCGAbiolinks:::getProjectSummary("TCGA-GBM")

# Search TCGA data
query <- GDCquery(
  project = c("TCGA-GBM"),
  data.category = c("Simple Nucleotide Variation"),
  access = c("open"))

# Download TCGA data 
#GDCdownload(query)

# Prepare data in SummarizedExperiment (default data structure for TCGAbiolinks) format. It can be saved as an R object.
## Mutation data is MAF files. I think these can not be saved as SummarizedExperiment.

#data = GDCprepare(query,
#                  summarizedExperiment = TRUE,
#                  save = TRUE,
#                  save.filename = gsub(":", "_", paste0(query$project, date(), "GDCprepare.rda"))) 


# Load GDC prepare file as 'data'
load("TCGA-GBMWed Oct 16 12_14_57 2024GDCprepare.Rdata")

# Check metadata
#TODO: Problem with receiving the data in summarizedExperiment format. Without it, clinical data is not accessible with colData function.
#colData(data)

barcode_list_for_metadata = unique(data$Tumor_Sample_Barcode)
metadata = as.data.frame(colDataPrepare(barcode_list_for_metadata))
metadata$barcode = gsub("-", ".", metadata$barcode) # using hyphens or dash sign (in general special characters) are not allowed in variable names. 
rownames(metadata) = metadata$barcode
# Write metadata
#saveRDS(metadata, file = gsub(":", "_", paste0("metadata/",query$project, date(), "_metadata.rds")))

# Read maf file
maf = data %>% maftools::read.maf()

## Maf Summary Visualization
data.frame(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

# Convert maf to binary matrix
maf_binary = maf2matrix(data)

# save as tsv
#fwrite(maf_binary, gsub(":", "_", paste0("binarymatrix/",query$project, date(), "_binary_mat.tsv")), sep = "\t", row.names = TRUE)

