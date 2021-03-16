##############################################################################
## Title: Robustness Speliote SNP Matching
## Author: Regina Manansala
## Date Created: 9-September-2019
## Date Modified: 9-September-2019
##############################################################################

library(data.table)
library(dplyr)

speliotes <- fread("~/Downloads/speliotes_bmi_imp.txt")
speliotes_bmi <- fread("~/Downloads/Speliotes_2010_BMI_SNPs.txt", header = F)

speliotes_nomatch <- speliotes_bmi[!(speliotes_bmi$V1 %in% spel_snps),]

write.table(speliotes_nomatch, "~/Documents/UKBB Income/BMI_Geno_Files/Speliotes_NoMatch.txt", sep = "\t", row.names = F, col.names = F, quote = F)
