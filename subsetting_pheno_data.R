library(dplyr)
library(tidyverse)

chr20genes <- read_table("chr20_GD462.signalGeneQuantRPKM_plink.txt")
targetgenes <- unlist(as.vector(read.table("chr10genes.txt")))

included <- colnames(chr10genes)[colnames(chr10genes) %in% targetgenes]

ID <- chr10genes[,1:2]
Genes_to_include <- chr10genes[,included]
subsetted_genes <- cbind(ID, Genes_to_include)

cat("IS NA IN MATRIX?")
TRUE %in% is.na(subsetted_genes)

write.table(subsetted_genes, "subsetted_genes.txt", quote = F, sep = " ", row.names = F)
