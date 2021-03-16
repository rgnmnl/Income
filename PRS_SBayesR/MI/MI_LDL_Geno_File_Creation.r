##############################################################################
## Title: MI LDL Genotype File Creation
## Version: 1
## Note: For SAIGE gwas results (used all 150k samples)
## Author: Regina Manansala
## Date Created: 15-March-2021
## Date Modified: 15-March-2021
##############################################################################

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)
library(tibble)

## Read LDL SNPs
# all_snps <- fread("/raid-04/SPH/pauer/manansa2/inc_prs/all_snps.txt")
ldl_snps <- fread("/raid-04/SPH/pauer/manansa2/inc_prs/SNP_List/MI_LDL_SNP_list.txt", header = FALSE)

## Add LDL snps to SNP reference file
snps <- fread("SNP_List/SNP_reference.txt")
names(ldl_snps) <- names(snps)
snps <- rbind(snps, ldl_snps)
write.table(snps, "SNP_List/SNP_reference.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

## Read BED files
geno_path <- "/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/LD_BED"
chr1 <- BEDMatrix(paste0(geno_path, "/ldsub_chr1.bed"))
chr2 <- BEDMatrix(paste0(geno_path, "/ldsub_chr2.bed"))
chr3 <- BEDMatrix(paste0(geno_path, "/ldsub_chr3.bed"))
chr4 <- BEDMatrix(paste0(geno_path, "/ldsub_chr4.bed"))
chr5 <- BEDMatrix(paste0(geno_path, "/ldsub_chr5.bed"))
chr6 <- BEDMatrix(paste0(geno_path, "/ldsub_chr6.bed"))
chr7 <- BEDMatrix(paste0(geno_path, "/ldsub_chr7.bed"))
chr8 <- BEDMatrix(paste0(geno_path, "/ldsub_chr8.bed"))
chr9 <- BEDMatrix(paste0(geno_path, "/ldsub_chr9.bed"))
chr10 <- BEDMatrix(paste0(geno_path, "/ldsub_chr10.bed"))
chr11 <- BEDMatrix(paste0(geno_path, "/ldsub_chr11.bed"))
chr12 <- BEDMatrix(paste0(geno_path, "/ldsub_chr12.bed"))
chr13 <- BEDMatrix(paste0(geno_path, "/ldsub_chr13.bed"))
chr14 <- BEDMatrix(paste0(geno_path, "/ldsub_chr14.bed"))
chr15 <- BEDMatrix(paste0(geno_path, "/ldsub_chr15.bed"))
chr16 <- BEDMatrix(paste0(geno_path, "/ldsub_chr16.bed"))
chr17 <- BEDMatrix(paste0(geno_path, "/ldsub_chr17.bed"))
chr18 <- BEDMatrix(paste0(geno_path, "/ldsub_chr18.bed"))
chr19 <- BEDMatrix(paste0(geno_path, "/ldsub_chr19.bed"))
chr20 <- BEDMatrix(paste0(geno_path, "/ldsub_chr20.bed"))
chr21 <- BEDMatrix(paste0(geno_path, "/ldsub_chr21.bed"))

## Subset BED files to LDL SNPs only
geno <- wg[, colnames(wg) %in% c(ldl_snps$V2, ldl_snps$V3)]
wg_sub <- geno %>% as.data.frame() %>% rownames_to_column('ID') %>% mutate(ID = gsub("_.*","",ID) %>% as.integer())

## Read phenotype file 
ukbb_mi <- fread("ukbb_prs_pheno_w_MI.txt")

## Write to file
MI_SNP_Genotypes <- left_join(ukbb_mi[, c("ID", "MI", "LDL")], wg_sub, by = "ID")
write.table(MI_SNP_Genotypes, "BOOTSTRAP/MI/MI_SNP_Genotypes.txt", sep = "\t", row.names = FALSE, col.names= TRUE, quote = FALSE)

############################################################
#################### Get MI & LDL BETAS ####################
############################################################

## Get build and test subsets
build <- fread("build_covar_gcta.txt")
mi_build <- MI_SNP_Genotypes[MI_SNP_Genotypes$ID %in% build$V1,]
mi_test <- MI_SNP_Genotypes[!(MI_SNP_Genotypes$ID %in% build$V1),]

## Run regression and save betas and se
build_betas <- data.frame(matrix(ncol=4, nrow=0))
names(build_betas) <- c("build_mi_beta", "build_mi_se", "build_ldl_beta", "build_ldl_se")
# names(build_betas) <- c("test_mi_beta", "test_mi_se", "test_ldl_beta", "test_ldl_se", "build_mi_beta", "build_mi_se", "build_ldl_beta", "build_ldl_se")

for(i in names(mi_build)[4:30]){
# 	build_betas[i, "test_mi_beta"] <- summary(lm(MI ~ get(i), data=mi_test ))$coef[2,1]
# 	build_betas[i, "test_mi_se"] <- summary(lm(MI ~ get(i), data=mi_test ))$coef[2,2]
# 	build_betas[i, "test_ldl_beta"] <- summary(lm(LDL ~ get(i), data=mi_test ))$coef[2,1]
# 	build_betas[i, "test_ldl_se"] <- summary(lm(LDL ~ get(i), data=mi_test ))$coef[2,2]
	build_betas[i, "build_mi_beta"] <- summary(lm(MI ~ get(i), data=mi_build ))$coef[2,1]
	build_betas[i, "build_mi_se"] <- summary(lm(MI ~ get(i), data=mi_build ))$coef[2,2]
	build_betas[i, "build_ldl_beta"] <- summary(lm(LDL ~ get(i), data=mi_build ))$coef[2,1]
	build_betas[i, "build_ldl_se"] <- summary(lm(LDL ~ get(i), data=mi_build ))$coef[2,2]
}

build_betas <- build_betas %>% rownames_to_column("SNP") %>% mutate(SNP = gsub("_.*","",SNP))

## Write to file
MI_LDL_SNPs <- fread("BOOTSTRAP/MI/MI_LDL_SNPs.txt")
MI_LDL_SNPs <- left_join(MI_LDL_SNPs, build_betas, by = "SNP")
write.table(MI_LDL_SNPs[, -c("SNP_1", "SNP_2")], "BOOTSTRAP/MI/MI_LDL_SNPs.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

############################################################
######################## Get MI PRS ########################
############################################################

## Get MI PRS files
mi_files <- list.files(path="/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n/score", pattern = ".profile")

## Import CHR PRS files
for(i in mi_files){
	temp <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n/score/",i)) #%>% right_join(., sample.ids, by = c("FID" = "V1"))
	assign(gsub("\\..*","", i), temp)
}

## Combine CHR PRS to a single file
mi_full <- mi_chr1

for(k in c(2:22)){
chri <- get(paste0("mi_chr", k))
	names(chri)[6] <- paste0("SCORE_", k)
	mi_full <- left_join(mi_full, chri[, c(2,6)], by = "IID")
}
head(mi_full)

names(mi_full)[6] <- "SCORE_1"

## Write to file
MI_PRS <- left_join(ukbb_mi[, c("ID", "MI", "LDL")], mi_full[, c(2, 6:27)], by = c("ID" = "IID"))
write.table(MI_PRS, "BOOTSTRAP/MI/MI_PRS.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

#############################################################
######################## Get LDL PRS ########################
#############################################################

## Get LDL PRS files
ldl_files <- list.files(path="/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/LDL/LDL_qc_n/score", pattern = ".profile")

## Import CHR PRS files
for(i in ldl_files){
	temp <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/LDL/LDL_qc_n/score/",i)) #%>% right_join(., sample.ids, by = c("FID" = "V1"))
	assign(gsub("\\..*","", i), temp)
}

## Combine CHR PRS to a single file
ldl_full <- ldl_chr1

for(k in c(2:22)){
chri <- get(paste0("ldl_chr", k))
	names(chri)[6] <- paste0("SCORE_", k)
	ldl_full <- left_join(ldl_full, chri[, c(2,6)], by = "IID")
}
head(ldl_full)

names(ldl_full)[6] <- "SCORE_1"

## Write to file
LDL_PRS <- left_join(ukbb_mi[, c("ID", "MI", "LDL")], ldl_full[, c(2, 6:27)], by = c("ID" = "IID"))
write.table(LDL_PRS, "BOOTSTRAP/LDL/LDL_PRS.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
















