##############################################################################
## Title: SEM Robustness BED file
## Author: Regina Manansala
## Date Created: 10-September-2019
## Date Modified: 16-September-2019
##############################################################################

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)

#######################

geno_path <- "~/Data/UKBB_PH/BED_BMI_Robustness/"

chr1 <- BEDMatrix(paste0(geno_path, "chr1_BMI_rbst.bed"))
chr2 <- BEDMatrix(paste0(geno_path, "chr2_BMI_rbst.bed"))
chr3 <- BEDMatrix(paste0(geno_path, "chr3_BMI_rbst.bed"))
chr4 <- BEDMatrix(paste0(geno_path, "chr4_BMI_rbst.bed"))
chr5 <- BEDMatrix(paste0(geno_path, "chr5_BMI_rbst.bed"))
chr6 <- BEDMatrix(paste0(geno_path, "chr6_BMI_rbst.bed"))
chr7 <- BEDMatrix(paste0(geno_path, "chr7_BMI_rbst.bed"))
chr8 <- BEDMatrix(paste0(geno_path, "chr8_BMI_rbst.bed"))
chr9 <- BEDMatrix(paste0(geno_path, "chr9_BMI_rbst.bed"))
chr10 <- BEDMatrix(paste0(geno_path, "chr10_BMI_rbst.bed"))
chr11 <- BEDMatrix(paste0(geno_path, "chr11_BMI_rbst.bed"))
chr12 <- BEDMatrix(paste0(geno_path, "chr12_BMI_rbst.bed"))
chr13 <- BEDMatrix(paste0(geno_path, "chr13_BMI_rbst.bed"))
chr14 <- BEDMatrix(paste0(geno_path, "chr14_BMI_rbst.bed"))
chr15 <- BEDMatrix(paste0(geno_path, "chr15_BMI_rbst.bed"))
chr16 <- BEDMatrix(paste0(geno_path, "chr16_BMI_rbst.bed"))
chr17 <- BEDMatrix(paste0(geno_path, "chr17_BMI_rbst.bed"))
chr18 <- BEDMatrix(paste0(geno_path, "chr18_BMI_rbst.bed"))
chr19 <- BEDMatrix(paste0(geno_path, "chr19_BMI_rbst.bed"))
chr20 <- BEDMatrix(paste0(geno_path, "chr20_BMI_rbst.bed"))
chr21 <- BEDMatrix(paste0(geno_path, "chr21_BMI_rbst.bed"))

wg <- ColumnLinkedMatrix(chr1, chr2, chr3, chr4, chr5, chr6, chr7,
						chr8, chr9, chr10, chr11, chr12, chr13, chr14,
						chr15, chr16, chr17, chr18, chr19, chr20, chr21)
geno <- as.matrix(wg) %>% as.data.frame()
geno$ID <- gsub("_.*","",rownames(geno)) %>% as.integer()
foo <- gsub("_.*", "", colnames(geno))
foo[43] <- paste0(foo[43], ".1")
names(geno) <- foo

########################

geno_path <- "~/Data/UKBB_PH/INC/SEM_FILES/GENO_FILES/"

chr1 <- BEDMatrix(paste0(geno_path, "sem_chr1.bed"))
chr2 <- BEDMatrix(paste0(geno_path, "sem_chr2.bed"))
chr3 <- BEDMatrix(paste0(geno_path, "sem_chr3.bed"))
chr4 <- BEDMatrix(paste0(geno_path, "sem_chr4.bed"))
chr5 <- BEDMatrix(paste0(geno_path, "sem_chr5.bed"))
chr6 <- BEDMatrix(paste0(geno_path, "sem_chr6.bed"))
chr7 <- BEDMatrix(paste0(geno_path, "sem_chr7.bed"))
chr8 <- BEDMatrix(paste0(geno_path, "sem_chr8.bed"))
chr9 <- BEDMatrix(paste0(geno_path, "sem_chr9.bed"))
chr10 <- BEDMatrix(paste0(geno_path, "sem_chr10.bed"))
chr11 <- BEDMatrix(paste0(geno_path, "sem_chr11.bed"))
chr12 <- BEDMatrix(paste0(geno_path, "sem_chr12.bed"))
chr13 <- BEDMatrix(paste0(geno_path, "sem_chr13.bed"))
chr14 <- BEDMatrix(paste0(geno_path, "sem_chr14.bed"))
chr15 <- BEDMatrix(paste0(geno_path, "sem_chr15.bed"))
chr16 <- BEDMatrix(paste0(geno_path, "sem_chr16.bed"))
chr17 <- BEDMatrix(paste0(geno_path, "sem_chr17.bed"))
chr18 <- BEDMatrix(paste0(geno_path, "sem_chr18.bed"))
chr19 <- BEDMatrix(paste0(geno_path, "sem_chr19.bed"))
chr20 <- BEDMatrix(paste0(geno_path, "sem_chr20.bed"))
chr21 <- BEDMatrix(paste0(geno_path, "sem_chr21.bed"))

wg <- ColumnLinkedMatrix(chr1, chr2, chr3, chr4, chr5, chr6, chr7,
						chr8, chr9, chr10, chr11, chr12, chr13, chr14,
						chr15, chr16, chr17, chr18, chr19, chr20, chr21)
geno_sem <- as.matrix(wg) %>% as.data.frame()
geno_sem$ID <- gsub("_.*","",rownames(geno)) %>% as.integer()

##################

fam <- "sem_all_chrom_rn.fam"
bed <- "sem_all_chrom.bed"
bim <- "sem_all_chrom.bim"

plink <- read.plink(bed, bim, fam)

geno_plink <- as(plink$genotypes, "numeric")
geno_plink$ID <- plink$fam$member
geno_plink <- as.data.frame(geno_plink)

#################

geno_sem_sub <- geno_sem[order(geno_sem$ID),colnames(geno_sem) %in% colnames(geno)]
geno_sub <- geno[order(geno$ID), colnames(geno) %in% colnames(geno_sem_sub)]
geno_plink_sub <- geno_plink[order(geno_plink$ID), colnames(geno_plink) %in% colnames(geno_sem_sub)]


head(geno_plink_sub[1:10, 1:10])
head(geno_sem_sub[1:10, 1:10])
head(geno_sub[1:10, 1:10])





