
library(data.table)
library(dplyr)

#oad("D:/GitHub/sarok01/tcga_data/GDCprepare/TCGACOAD_20241030_GDCprepare.Rdata")
load("D:/GitHub/sarok01/tcga_data/GDCprepare/TCGAREAD_20241030_GDCprepare.Rdata")

#unique_mutations = data %>% distinct(Chromosome, Start_Position, End_Position)

# Genes of interest from Fisher results
#subset_gene_list = c("KRAS","TP53","PIK3CA","DCST1","TTLL4","ST3GAL6","GOT1L1","NOLC1","ALDH3A1", "LDLRAP1" ,"FAM20B" , "RBPMS",
#                                    "GPT","MAN1B1","VPS37C","NOX5","POFUT2","SEPTIN5","SOAT1","SCAP","AP1S1","JADE3","VPS39","PPP1R3F",
#                                    "PDIA5","ARID4B","IRAK4","FAM209A","XK","STK11","LYAR","GPR22","GAS2L2","RPP40",
#                                    "VPS35L","PMFBP1","SLC1A7","SAFB","CXADR","UBQLN2","AADACL3","USP42","CD180","KLHL14","GPC4","FCHSD2","BICRA","APC")

# Keep only interesting genes for faster computation
#data_subset = subset(data, Hugo_Symbol %in% subset_gene_list)
data_subset = data

# Get high and moderate impact mutations
data_subset_highimpact = subset(data_subset, IMPACT == "HIGH" | IMPACT == "MODERATE")

# calculate extract the numeric information from polyphen and sift columns
data_subset_highimpact$polyphen_numeric = lapply(data_subset_highimpact$PolyPhen, function(x) as.numeric(strsplit(as.character(x), split = "\\(|\\)")[[1]][2]))
data_subset_highimpact$sift_numeric = lapply(data_subset_highimpact$SIFT, function(x) as.numeric(strsplit(as.character(x), split = "\\(|\\)")[[1]][2]))

# Only Polyphen subsetting
data_subset_highimpact_polyphen = subset(data_subset_highimpact, is.na(polyphen_numeric) | polyphen_numeric > 0.446)

# Only Sift subsetting
data_subset_highimpact_sift = subset(data_subset_highimpact, is.na(sift_numeric) | sift_numeric < 0.05)

# SIFT and PolyPhen Filtering
data_subset_highimpact_polyphenSIFT = subset(data_subset_highimpact, (is.na(sift_numeric) | sift_numeric < 0.05) & (is.na(polyphen_numeric) | polyphen_numeric > 0.446))


# Save 
saveRDS(data_subset_highimpact_polyphenSIFT, "TCGAREAD_20241030_GDCprepare_PolyphenSIFTfiltered.Rdata")
# Access via API - multiple entries

# library(httr)
# library(jsonlite)
# library(xml2)

# server <- "https://rest.ensembl.org"
# ext <- "/vep/human/hgvs"
# r <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = '{ "hgvs_notations" : ["ENST00000378156:c.4024G>A", "ENST00000370079:c.940C>T"] }')

# stop_for_status(r)

# # use this if you get a simple nested list back, otherwise inspect its structure
# # head(data.frame(t(sapply(content(r),c))))
# head(fromJSON(toJSON(content(r))))