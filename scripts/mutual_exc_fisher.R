library(data.table) # fread
library(ggplot2)

#maf_binary2 = data.frame(fread("Desktop/r2/mutual_exclusivity/R_approach/gbm_binary.tsv")) #row names are on the first column since we converted from MAF to TSV

#row.names(maf_binary2) = maf_binary2$V1

#maf_binary2 = maf_binary2[,-1]

## duplicated columns
#asd = colnames(maf_binary2)[duplicated(colnames(maf_binary2))]

#colnames(maf_binary2)[which(colnames(maf_binary2) == "TCGA.27.2526.01")]

#duplicated_col_ex = maf_binary2[,26:27]

#View(duplicated_col_ex[which(duplicated_col_ex[[1]] != duplicated_col_ex[[2]]),]) # discrepancy in 3 genes in this example


# remove duplicates because it affects filtering (due to over representation in the row sum)
setwd("D:/R_projects/mutual_exc/tcga_data/")
maf_binary = fread("binarymatrix/TCGA-GBMWed Oct 16 14_29_56 2024_binary_mat.tsv")

notdup_col_pos = which(!duplicated(colnames(maf_binary)))

maf_binary = data.frame(maf_binary)[,notdup_col_pos]

# add row names

row.names(maf_binary) = maf_binary$V1

maf_binary = maf_binary[,-1]

### filtering

## row sum
maf_binary$row.sum = rowSums(maf_binary)

## plot row sum

ggplot(maf_binary, aes(x = row.sum)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black")

## remove rowsum < 5 

maf_binary_filtered = maf_binary[which(maf_binary$row.sum >= 5),]

### Sampling

# Calculate mutation frequency for each gene
mutation_frequencies <- maf_binary_filtered$row.sum

# Normalize frequencies to get proportions
normalized_frequencies <- mutation_frequencies / sum(mutation_frequencies)
print(normalized_frequencies)

# Function to sample weighted gene pairs
weighted_gene_pairs <- function(genes, weights, num_pairs) {
  sampled_pairs <- list()
  
  for (i in 1:num_pairs) {
    # Sample two genes with weights
    sampled_genes <- sample(genes, 2, prob = weights, replace = FALSE)
    sampled_pairs[[i]] <- sampled_genes
  }
  
  return(sampled_pairs)
}

# List of genes
genes <- rownames(maf_binary_filtered)

# Number of gene pairs to sample
num_pairs <- 10

# Perform weighted sampling
sampled_gene_pairs <- weighted_gene_pairs(genes, normalized_frequencies, num_pairs)
print(sampled_gene_pairs)

### Pair-wise testing

maf_binary_t = t(maf_binary_filtered[,-length(colnames(maf_binary_filtered))]) #transposed and removed rowsum 

# Function to calculate Fisher's Exact Test for mutual exclusivity

look_up_table = list()
#look_up_table[["290-76-9-7"]] = 0.055363

mutual_exclusivity_test <- function(gene1, gene2, transposed_df, look_up_table) {
  # Create a contingency table
  contingency_table <- table(transposed_df[,gene1], transposed_df[,gene2])
  print(contingency_table)
  
  table_values_str = paste(contingency_table, collapse = "-") # first number x axis - second number y axis. number order: t00 - t01 - t10 - t11
  
  # Apply Fisher's exact test
  if (!is.null(look_up_table[[table_values_str]])) { # if key exists, get p-value from the there
    
    return(list(p_value = look_up_table[[table_values_str]], contingency_table = contingency_table, look_up_table = look_up_table))
    
  } else { # if key does not exist, run fisher test 
    fisher_test_result <- fisher.test(contingency_table)
    look_up_table[[table_values_str]] = fisher_test_result$p.value
    return(list(p_value = fisher_test_result$p.value, contingency_table = contingency_table, look_up_table = look_up_table))
  }
}

# Example: Testing for mutual exclusivity between Gene 1 and Gene 2

gene1 = "PTCH1"
gene2 = "ZNF81"

maf_binary_t[,gene2][maf_binary_t[,gene2] == 1] <- 2
result <- mutual_exclusivity_test(gene1, gene2, maf_binary_t, look_up_table)

look_up_table = result$look_up_table

# View the result
#print(result$p_value)
#print(result$contingency_table)


# Loop through all pairs of genes
#res_pval = list()
  
#for (i in 1:length(sampled_gene_pairs)) {
#    gene1 = sampled_gene_pairs[[i]][1]
#    gene2 = sampled_gene_pairs[[i]][2]
    
#    result <- mutual_exclusivity_test(gene1, gene2, maf_binary_t, look_up_table)
#    look_up_table = result$look_up_table
#    res_pval[[paste(gene1,gene2, collapse = "-")]] = result$p_value
#    cat("Mutual exclusivity between", gene1, "and", gene2, ": p-value =", result$p_value, "\n")
#}



### 3 genes

sep_gene = "TP53" #seperator gene

mutated_matrix = maf_binary_t[which(maf_binary_t[,sep_gene] == 1),]

unmutated_matrix = maf_binary_t[which(maf_binary_t[,sep_gene] == 0),]


res_pval_3gene = data.frame(
                  "gene_pair" = character(),
                  "mutated_pval" = numeric(),
                  "unmutated_pval" = numeric(),
                  "mut_contable" = character(),
                  "unmut_contable" = character()
)

#for (i in 1:length(sampled_gene_pairs)) {
#  gene1 = sampled_gene_pairs[[i]][1]
#  gene2 = sampled_gene_pairs[[i]][2]
  
#  # fisher on mutated
#  result_mutated <- mutual_exclusivity_test(gene1, gene2, mutated_matrix, look_up_table)
#  look_up_table = result_mutated$look_up_table
#  # fisher on unmutated
#  result_unmutated <- mutual_exclusivity_test(gene1, gene2, unmutated_matrix, look_up_table)
#  look_up_table = result_unmutated$look_up_table
  
#  # append p values
#  pval_3gene_calculated = data.frame(paste(gene1,gene2, collapse = "-"), result_mutated$p_value, result_unmutated$p_value)
#  colnames(pval_3gene_calculated) = c("gene_pair", "mutated_pval", "unmutated_pval")
#  res_pval_3gene = rbind(res_pval_3gene, pval_3gene_calculated)  # gene pair name, mutated pval, unmutated pval
#  cat("Mutual exclusivity between", gene1, "and", gene2, ": mutated p-value =", result_mutated$p_value, "unmutated p-value =", result_unmutated$p_value,"\n")
#} 

iteration_counter = 0
for (i in 1:(ncol(mutated_matrix) - 1)) { #mutated and unmutated matrices have the ncol since they have the same genes.
  for (j in (i + 1):ncol(mutated_matrix)) {
    gene1 = colnames(mutated_matrix)[i]
    gene2 = colnames(mutated_matrix)[j]
    
    iteration_counter = iteration_counter + 1
    if (iteration_counter %% 1000 == 0) {
      print(paste0("Iteration ", iteration_counter))
    }
    
    if (sum(mutated_matrix[,i]) %in% c(0, nrow(mutated_matrix)) || 
        sum(unmutated_matrix[,i]) %in% c(0, nrow(unmutated_matrix)) ||
        sum(mutated_matrix[,j]) %in% c(0, nrow(mutated_matrix)) || 
        sum(unmutated_matrix[,j]) %in% c(0, nrow(unmutated_matrix))) {
      next # skip iteration if a gene has only 0 or 1 values - can not create a 2x2 contingency table
    } else {
    # fisher on mutated
    result_mutated <- mutual_exclusivity_test(gene1, gene2, mutated_matrix, look_up_table)
    look_up_table = result_mutated$look_up_table
    # fisher on unmutated
    result_unmutated <- mutual_exclusivity_test(gene1, gene2, unmutated_matrix, look_up_table)
    look_up_table = result_unmutated$look_up_table

    # append p values
    pval_3gene_calculated = data.frame(paste(gene1,gene2, collapse = "-"),
                                       result_mutated$p_value, 
                                       result_unmutated$p_value, 
                                       paste(result_mutated$contingency_table, collapse = "-"),
                                       paste(result_unmutated$contingency_table, collapse = "-"))
    colnames(pval_3gene_calculated) = c("gene_pair", "mutated_pval", "unmutated_pval", "mut_contable", "unmut_contable")
    res_pval_3gene = rbind(res_pval_3gene, pval_3gene_calculated)  # gene pair name, mutated pval, unmutated pval
    #cat("Mutual exclusivity between", gene1, "and", gene2, ": mutated p-value =", result_mutated$p_value, "unmutated p-value =", result_unmutated$p_value,"\n")
    }
  }
}

