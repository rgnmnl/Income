##############################################################################
## Title: UKBB Income PLINK GWAS Results Analysis
## Author: Regina Manansala
## Date Created: 08-January-2019
## Date Modified: 11-February-2019
##############################################################################

library(data.table)
library(dplyr)
library(qqman)
library(manhattanly)
library(tidyr)
library(stringr)

## Read plink results
inc_regression_files <- list.files(pattern="inc_chr(\\d|\\d\\d)_assoc.Income.glm.logistic$")
inc_regression_files <- list.files(pattern="inc_chr(\\d|\\d\\d)_assoc.Income.glm.linear$")
inc_results <- do.call(rbind, lapply(inc_regression_files, function(x) fread(x, stringsAsFactors = FALSE)))

## MAF filter
inc_maf_filt <- inc_results[inc_results$A1_FREQ >= .01,]

## Read lists of SNPs
bmi_snp <- fread("../BMI_SNPs.txt")
t2d_snp <- fread("../T2D_SNPs.csv")
cad_snp <- fread("../CAD_SNPs.csv")

# bmi_snp <- read.table("../bmi_grch37_ordered.txt", header=T)

## Get GWAS results for specific SNPs
inc_bmi <- subset(inc_results, ID %in% bmi_snp$SNP)
# inc_bmi <- subset(inc_results, POS %in% bmi_snp$POS)
inc_t2d <- subset(inc_results, ID %in% t2d_snp$SNP)
inc_cad <- subset(inc_results, ID %in% cad_snp$SNP)

## Export
write.table(inc_bmi, "inc_bmi.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(inc_t2d, "inc_t2d.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(inc_cad, "inc_cad.txt", sep="\t", row.names=F, col.names=T, quote=F)

write.table(inc_bmi, "con_inc_bmi.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(inc_t2d, "con_inc_t2d.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(inc_cad, "con_inc_cad.txt", sep="\t", row.names=F, col.names=T, quote=F)

## Plot results
names(inc_maf_filt)[1] <- "CHROM"
plot <- manhattanr(inc_maf_filt, chr="CHROM", snp="ID", bp="POS", p="P")

jpeg("manhattan.jpg", height=5, width=7, units='in', quality=100, res=500)
manhattan(inc_maf_filt, chr="CHROM", snp="ID", bp="POS", p="P")
dev.off()

jpeg("qq.jpg", height=5, width=5, units='in', quality=100, res=500)
qq(inc_maf_filt$P, main="QQ Plot")
dev.off()

bmi_alleles <- str_split_fixed(bmi_snp$Alleles, '/', 2)
bmi_snp$REF <- bmi_alleles[,1]
bmi_snp$ALT1 <- bmi_alleles[,2]

write.table(inc_bmi, "inc_bmi_plink2.txt", sep="\t", row.names=F, col.names=T, quote=F)
