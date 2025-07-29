# Load packages
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(SummarizedExperiment)

TCGAbiolinks:::getProjectSummary("TCGA-COAD")
TCGAbiolinks:::getProjectSummary("TCGA-READ")

# Set directory
setwd("D:/GitHub/sarok01/tcga_data/")

# Function for getting the TCGA gene expression data. Th
get_gene_exp = function(project_name = c("TCGA-COAD"),
                        data_cat = c("Transcriptome Profiling"), 
                        type = "gene_expression", # under transcriptome profiling, there are also miRNA samples 
                        sample_type = "Primary Tumor") #  there might be different types of samples: primary tumors, normal tissue ...
                                                          {
  query = GDCquery(
                   project = project_name,
                   data.category = data_cat,
                   access = c("open"))

  query[[1]][[1]] = query[[1]][[1]][which(query[[1]][[1]][["type"]] == 'gene_expression' & query[[1]][[1]][["sample_type"]] == 'Primary Tumor'),]

  # Download TCGA data 
  GDCdownload(query)

  data = GDCprepare(query,
                    summarizedExperiment = TRUE,
                    save = TRUE,
                    save.filename = gsub("-", "", paste0(paste(project_name, collapse = ""),'_',Sys.Date(), "_GDCprepare.Rdata"))) 
  return(data)
}

# Get RNAseq data for COAD and READ
data = get_gene_exp(project_name = c("TCGA-COAD", "TCGA-READ"), data_cat = c("Transcriptome Profiling"))

rna_samples_coad = colnames(assay(data))

rna_samples_coad = gsub("-", ".", rna_samples_coad)

# Compare labels to Arezo's

cms_pred = read.csv(file = "D:/GitHub/sarok01/cms_subtyping/data/cms_prediction_nanocmser_arezo.txt", sep = " ")

length(intersect(rna_samples_coad, rownames(cms_pred)))

# Compare labels to mut exc binary matrix
binary = read.csv("binarymatrix/TCGACOADREAD_20241030_PolyPhenSIFTfiltered20241211_binarymat.tsv", sep = "\t")

length(intersect(rna_samples_coad, colnames(binary)))

# Function to extract first 3 dot-separated parts
get_first_three <- function(x) {
  parts <- strsplit(x, "\\.")[[1]]
  paste(parts[1:3], collapse = ".")
}

get_first_four <- function(x) {
  parts <- strsplit(x, "\\.")[[1]]
  paste(parts[1:4], collapse = ".")
}

get_first =  function(x) {
  parts <- strsplit(x, "\\.")[[1]]
  parts[1]
}

# Analyze participants/samples
rna_samples_coad_participant = sapply(rna_samples_coad, get_first_three)

rna_samples_coad_samplevial = sapply(rna_samples_coad, get_first_four)

rna_samples_id_df = data.frame("full_barcode" = rna_samples_coad, "first_three" = rna_samples_coad_participant, "first_four" = rna_samples_coad_samplevial)


binary_labels_samplevial = sapply(colnames(binary)[2:616], get_first_four)

length(intersect(rna_samples_coad_samplevial, binary_labels_samplevial))

length(intersect(sapply(rownames(cms_pred), get_first_four), binary_labels_samplevial))

binary_labels_participant = sapply(colnames(binary)[2:616], get_first_three)

length(intersect(rna_samples_coad_participant, binary_labels_participant))
length(intersect(sapply(rownames(cms_pred), get_first_three), binary_labels_participant))

# TODO:Save summarizedExperiment object.

# Ensembl ID remove version 

count_matrix = assay(data)

ensembl_ids_noversion = sapply(rownames(count_matrix), get_first)
names(ensembl_ids_noversion) = NULL

length(unique(ensembl_ids_noversion))

length(rownames(count_matrix)) == length(ensembl_ids_noversion)

rownames(count_matrix) = ensembl_ids_noversion






