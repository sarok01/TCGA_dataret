# Load packages
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(SummarizedExperiment)


# Set directory
setwd("D:/GitHub/sarok01/tcga_data/")
source(file = "scripts/maf2matrix.R")

# TCGA
TCGAbiolinks:::getProjectSummary("TCGA-READ")
TCGAbiolinks:::getProjectSummary("TCGA-COAD")

# Search TCGA data

# check for project types: https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
query <- GDCquery(
  project = c("TCGA-READ"),
  data.category = c("Simple Nucleotide Variation"),
  access = c("open"))

# Download TCGA data 
GDCdownload(query)

# Prepare data in SummarizedExperiment (default data structure for TCGAbiolinks) format. It can be saved as an R object.
## Mutation data is MAF files. I think these can not be saved as SummarizedExperiment.

data = GDCprepare(query,
                  summarizedExperiment = TRUE,
                  save = TRUE,
                  save.filename = gsub("-", "", paste0(query$project,'_',Sys.Date(), "_GDCprepare.Rdata"))) 


# Load GDC prepare file as 'data'

load('GDCprepare/TCGACOAD_20241030_GDCprepare.Rdata')

# Check metadata
#TODO: Problem with receiving the data in summarizedExperiment format. Without it, clinical data is not accessible with colData function.
#colData(data)

# Metadata preparation

barcode_list_for_metadata = unique(data$Tumor_Sample_Barcode)
metadata = as.data.frame(colDataPrepare(barcode_list_for_metadata))
metadata$barcode = gsub("-", ".", metadata$barcode) # using hyphens or dash sign (in general special characters) are not allowed in variable names. 
rownames(metadata) = metadata$barcode
# Write metadata
saveRDS(metadata, file = gsub("-", "", paste0(query$project,'_',Sys.Date(), "_metadata.rds")))

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
#fwrite(maf_binary, gsub("-", "", paste0("binarymatrix/",query$project, Sys.Date(), "_binarymat.tsv")), sep = "\t", row.names = TRUE)

