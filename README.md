# plink related rscripts


## Linux commands
- SNP association with genes
```
# One SNP vs All genes
plink --bfile [genotype-filename] --pheno [phenotype-filename] --all-pheno --snps [snp-name] --linear --out [directory/filename]
# Apply filter_pval_snp_genes.R with directory of result files
```
- Gene association with SNPs
```
# One Gene vs All SNPs
plink --bfile [genotype-filename] --pheno [phenotype-filename] --pheno-name [genename] --snps [snp-name] --linear --out [directory/filename]

# One Gene vs list of SNPs
plink --bfile [genotype-filename] --pheno [phenotype-filename] --pheno-name [genename] --extract [snp-list-file] --linear --out [directory/filename]

# Sorting and appply cutoff with P-value
awk '$9 < 0.05 {print $0}' ./Results/[results-filename] | sort -k9 -g > output_filename.txt
```
- For LD plot (Association with SNPs)
```
plink --bfile [genotype-filename] --extract [SNP-list-file] --recode HV --out [directory/filename]
```
- Make ped for boxplot
```
plink —bfile [filename] --snps [snp] --recode —pheno [chr] --pheno-name [gene] --out [filename]
```
