##############################################################################
## Title: UKBB Height Analysis
## Author: Regina Manansala
## Date Created: 10-February-2019
## Date Modified: 03-March-2019
##############################################################################

library(data.table)
library(dplyr)
library(qqman)
library(manhattanly)
library(tidyr)
library(stringr)

inc_dir <- "~/Data/UKBB_PH/INC/INC_CONTINUOUS_20190118/"
finc_dir <- "~/Data/UKBB_PH/INC/INC_CONTINUOUS_20190118/INC_FEMALE/"
minc_dir <- "~/Data/UKBB_PH/INC/INC_CONTINUOUS_20190118/INC_MALE/"

#IMPORT INCOME PLINK RESULTS
inc_regression_files <- list.files(path=inc_dir, pattern="inc_chr(\\d|\\d\\d)_assoc.Income.glm.linear$")
inc_results <- do.call(rbind, lapply(inc_regression_files, function(x) fread(paste0(inc_dir,x), stringsAsFactors = FALSE)))

#FEMALE
inc_regression_files <- list.files(path=finc_dir, pattern="inc_chr(\\d|\\d\\d)_assoc.Income.glm.linear$")
inc_results <- do.call(rbind, lapply(inc_regression_files, function(x) fread(paste0(finc_dir,x), stringsAsFactors = FALSE)))

#MALE
inc_regression_files <- list.files(path=minc_dir, pattern="inc_chr(\\d|\\d\\d)_assoc.Income.glm.linear$")
inc_results <- do.call(rbind, lapply(inc_regression_files, function(x) fread(paste0(minc_dir,x), stringsAsFactors = FALSE)))

# inc_regression_files <- list.files(path="../../../INC_CONTINUOUS_20190118", pattern="inc_chr(\\d|\\d\\d)_assoc.Income.glm.linear$")
# inc_results <- do.call(rbind, lapply(inc_regression_files, function(x) fread(paste0("../../../INC_CONTINUOUS_20190118/",x), stringsAsFactors = FALSE)))
## RENAME COLUMNS IN INC RESULTS
names(inc_results)[1] <- "CHROM"
colnames(inc_results) <- paste("inc", colnames(inc_results), sep="_")

inc_maf_filt <- inc_results[inc_results$A1_FREQ >= .01,]


#IMPORT HEIGHT PLINK RESULTS
plink_ht_files <- list.files(pattern="inc_chr(\\d|\\d\\d)_assoc.Height.glm.linear$")
ht_results <- do.call(rbind, lapply(plink_ht_files, function(x) fread(x, stringsAsFactors = FALSE)))
#FEMALE
ht_female_files <- list.files(pattern="inc_chr(\\d|\\d\\d)_assoc.Height.glm.linear$")
ht_results <- do.call(rbind, lapply(ht_female_files, function(x) fread(x, stringsAsFactors = FALSE)))
#MALE
ht_male_files <- list.files(pattern="inc_chr(\\d|\\d\\d)_assoc.Height.glm.linear$")
ht_results <- do.call(rbind, lapply(ht_male_files, function(x) fread(x, stringsAsFactors = FALSE)))
## RENAME COLUMNS IN HT RESULTS
names(ht_results)[1] <- "CHROM"
colnames(ht_results) <- paste("ht", colnames(ht_results), sep="_")

# 
ht_snps <- fread("../../Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt") %>% subset(P < 5*10^(-10))
# 
# ht_plink_subset <- subset(ht_results, ID %in% ht_snps$SNP)
# names(ht_plink_subset)[1] <- "CHROM"
# ht <- left_join(inc_ht, height, by=c("ID"="SNP"))

plink_ht <- left_join(ht_results[ht_ID %in% ht_snps$SNP, c(1:5, 12, 15)], inc_results[inc_ID %in% ht_snps$SNP, c("inc_ID", "inc_BETA", "inc_P")], by=c("ht_ID"="inc_ID"))
plink_fem_ht <- left_join(ht_results[ht_ID %in% ht_snps$SNP, c(1:5, 12, 15)], inc_results[inc_ID %in% ht_snps$SNP, c("inc_ID", "inc_BETA", "inc_P")], by=c("ht_ID"="inc_ID"))
plink_male_ht <- left_join(ht_results[ht_ID %in% ht_snps$SNP, c(1:5, 12, 15)], inc_results[inc_ID %in% ht_snps$SNP, c("inc_ID", "inc_BETA", "inc_P")], by=c("ht_ID"="inc_ID"))


# CHECK NUMBER OF MISMATCHED ALLELES
# dim(ht[ht$REF != ht$Tested_Allele & !is.na(ht$Tested_Allele),])

# CREATE TRANSFORMED BETA VARIABLE
# ht$BETA_2 <- ifelse(ht$REF == ht$Tested_Allele & !is.na(ht$Tested_Allele), ht$BETA.x, (-1)*ht$BETA.x)

# head(ht[ht$REF != ht$Tested_Allele & !is.na(ht$Tested_Allele), c("REF", "ALT1", "BETA.x", "BETA_2", "Tested_Allele", "Other_Allele")], 10)
# 
# dim(ht[ht$BETA.x == ht$BETA_2 & !is.na(ht$BETA_2),])

write.table(plink_ht, "plink_ht.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(plink_fem_ht, "plink_fem_ht.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(plink_male_ht, "plink_male_ht.txt", sep="\t", row.names=F, col.names=T, quote=F)