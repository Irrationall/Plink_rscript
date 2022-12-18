library(dplyr)
library(tidyverse)

chr_genes <- read_table("../../../Pheno_data/chr2_GD462.signalGeneQuantRPKM_plink.txt")
targetgenes <- unlist(as.vector(read.table("chr2targetgenes.txt")))

included <- colnames(chr_genes)[colnames(chr_genes) %in% targetgenes]

ID <- chr_genes[,1:2]
Genes_to_include <- chr_genes[,included]
subsetted_genes <- cbind(ID, Genes_to_include)

cat("IS NA IN MATRIX?")
TRUE %in% is.na(subsetted_genes)

write.table(subsetted_genes, "chr2_subsetted_genes.txt", quote = F, sep = " ", row.names = F)
