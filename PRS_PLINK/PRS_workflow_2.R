##############################################################################
## Title: PRS Workflow Step 2
## Purpose: Calculate PLINK score
## Author: Regina Manansala
## Date Created: 10-August-2019
## Date Modified: 29-May-2020
##############################################################################

# library(data.table)
# library(dplyr)
# 
# ht_path <- "/raid-05/SPH/pauer/UKBB_PH/INC/HEIGHT_ANALYSIS/PLINK/" %>% print()
# 
# plink_ht_files <- list.files(path=ht_path, pattern="inc_chr(\\d|\\d\\d)_assoc.Height.glm.linear$")
# ht <- do.call(rbind, lapply(plink_ht_files, function(x) fread(paste0(ht_path,x), stringsAsFactors = FALSE))) %>% select(., ID, ALT1, A1, BETA) %>% subset(., ID %in% ht_snps$SNP)
# ht$EA <- ifelse(ht$BETA < 0, ht$ALT1, ht$A1)
# ht$BETA2 <- ifelse(ht$BETA < 0, ht$BETA*(-1), ht$BETA)
# write.table(ht[,c("ID", "EA", "BETA2")], "~/Data/UKBB_PH/INC/SEM_FILES/score_ht.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

library(data.table)
library(dplyr)

#fastgwa_ht_files <- list.files(path="~/Data/manansa2/inc_prs/HEIGHT", pattern=".fastGWA")
#ht <- do.call(rbind, lapply(fastgwa_ht_files, function(x) fread(paste0("~/Data/manansa2/inc_prs/HEIGHT/",x), stringsAsFactors = FALSE))) %>% 
ht <- fread("Height_Index_SNPS.txt") %>% 
	select(., SNP, A1, A2, BETA)
ht$EA <- ifelse(ht$BETA < 0, ht$A2, ht$A1)
ht$BETA2 <- ifelse(ht$BETA < 0, ht$BETA*(-1), ht$BETA)
write.table(ht[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_ht.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

#fastgwa_bmi_files <- list.files(path="~/Data/manansa2/inc_prs/BMI", pattern=".fastGWA")
#bmi <- do.call(rbind, lapply(fastgwa_bmi_files, function(x) fread(paste0("BMI/", x), stringsAsFactors=FALSE))) %>% 
bmi <- fread("BMI_Index_SNPS.txt") %>% 
	select(., SNP, A1, A2, BETA)
bmi$EA <- ifelse(bmi$BETA < 0, bmi$A2, bmi$A1)
bmi$BETA2 <- ifelse(bmi$BETA < 0, bmi$BETA*(-1), bmi$BETA)
write.table(bmi[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_bmi.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

#fastgwa_inc_files <- list.files(path="~/Data/manansa2/inc_prs/INC", pattern=".fastGWA")
#inc <- do.call(rbind, lapply(fastgwa_inc_files, function(x) fread(paste0("INC/", x), stringsAsFactors=FALSE))) %>% 
inc <- fread("INC_Index_SNPS.txt") %>% 
	select(., SNP, A1, A2, BETA)
inc$EA <- ifelse(inc$BETA < 0, inc$A2, inc$A1)
inc$BETA2 <- ifelse(inc$BETA < 0, inc$BETA*(-1), inc$BETA)
write.table(inc[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_inc.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

edubin <- fread("EDUbin_Index_SNPS.txt") %>% 
	select(., SNP, A1, A2, BETA)
edubin$EA <- ifelse(edubin$BETA < 0, edubin$A2, edubin$A1)
edubin$BETA2 <- ifelse(edubin$BETA < 0, edubin$BETA*(-1), edubin$BETA)
write.table(edubin[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_edubin.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

edunum <- fread("EDUnum_Index_SNPS.txt") %>% 
	select(., SNP, A1, A2, BETA)
edunum$EA <- ifelse(edunum$BETA < 0, edunum$A2, edunum$A1)
edunum$BETA2 <- ifelse(edunum$BETA < 0, edunum$BETA*(-1), edunum$BETA)
write.table(edunum[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_edunum.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

alc <- fread("ALC_Index_SNPS.txt") %>% 
	select(., SNP, A1, A2, BETA)
alc$EA <- ifelse(alc$BETA < 0, alc$A2, alc$A1)
alc$BETA2 <- ifelse(alc$BETA < 0, alc$BETA*(-1), alc$BETA)
write.table(alc[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_alc.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

# dep <- fread("DEP_Index_SNPS.txt") %>% 
dep <- fread("DEP_Index_SNPS_v2.txt") %>% 
	select(., SNP, A1, A2, BETA)
dep$EA <- ifelse(dep$BETA < 0, dep$A2, dep$A1)
dep$BETA2 <- ifelse(dep$BETA < 0, dep$BETA*(-1), dep$BETA)
write.table(dep[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_dep_v2.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

smk <- fread("SMK_Index_SNPS.txt") %>% 
	select(., SNP, A1, A2, BETA)
smk$EA <- ifelse(smk$BETA < 0, smk$A2, smk$A1)
smk$BETA2 <- ifelse(smk$BETA < 0, smk$BETA*(-1), smk$BETA)
write.table(smk[,c("SNP", "EA", "BETA2")], "~/Data/manansa2/inc_prs/PRS/score_smk.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

###### PLINK SCORE CALCULATION

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=score_bmi

#printenv
#exit

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink2 \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --pfile /raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/ldsub_chr${SLURM_ARRAY_TASK_ID} \ ## USING UNMERGERD BED FILES
#     --bfile ~/Data/UKBB_PH/INC/SEM_FILES/GENO_FILES/sem_all_chrom \ ## USING MERGED BED FILES -- REMEMBER TO REMOVE SBATCH --array FLAG ABOVE
    --score ~/Data/manansa2/inc_prs/PRS/height_prs_snps.txt \
    --out ht_score


##### READ PLINK SCORE OUTPUTS

library(data.table)
library(dplyr)

# files <- list.files(pattern=".sscore")
# ht <- do.call(rbind, lapply(files, function(x) fread(x), stringsAsFactors = FALSE)))
# ht <- do.call(rbind, lapply(files, function(x) fread(x) %>% assign(), stringsAsFactors = FALSE)))

pheno <- fread("/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/GWAS/build_pheno.txt")
pheno_test <- fread("/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/GWAS/test_pheno.txt")

files <- list.files(pattern=".sscore")
files <- gsub("\\..*","", files)

for(i in files){
	foo <- fread(paste0(i, ".sscore")) %>% subset(., IID %in% pheno$IID)
# 	assign(gsub("\\..*","", i), foo)
	write.table(foo, paste0(i, ".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
}

for(i in files){
	foo_test <- fread(paste0(i, ".sscore")) %>% subset(., IID %in% pheno_test$IID)
# 	assign(gsub("\\..*","", i), foo_test)
	write.table(foo_test, paste0(i, "_test", ".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
}

# test <- list.files(pattern="_test.txt")
# test <- gsub("\\..*","", test)
# for(i in test){
# 	test <- fread(paste0(i, ".txt"))
# 	assign(gsub("\\..*","", i), test)
# # 	write.table(foo_test, paste0(i, "_test", ".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
# }
