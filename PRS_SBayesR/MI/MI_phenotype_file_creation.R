##############################################################################
## Title: UKBB MI phenotype using ICD-10
## Version: 2
## Author: Regina Manansala
## Date Created: 11-March-2021
## Date Modified: 13-March-2021
##############################################################################

library(data.table)
library(dplyr)

## Import phenotype files
ukb_prs <- fread("inc_prs/ukbb_prs_pheno_w_t2d.txt")
test_ids <- fread("inc_prs_gctb/gctb/test_subset_IDs.txt") 
build <- fread("inc_prs/build_covar_gcta.txt")

## Import UKBB data
source('ukb23544_v2.r')
inc <- bd
source("ukb43095.r")

## Get diagnosis variables and subset by diagnosis
ukb_mi <- select(bd, c("f.eid", "f.2010.0.0", "f.2090.0.0", "f.41202.0.0", "f.41202.0.1", "f.41204.0.0", "f.41204.0.1"))

ukb_mi$primary <- ifelse(grepl("I21|I22|I23|I241|I252", ukb_mi$f.41202.0.0), 1, 0)
ukb_mi$primary_v2 <- ifelse(grepl("I21|I22|I23|I241|I252", ukb_mi$f.41202.0.1), 1, 0)

ukb_mi$secondary <- ifelse(grepl("I21|I22|I23|I241|I252", ukb_mi$f.41204.0.0), 1, 0)
ukb_mi$secondary_v2 <- ifelse(grepl("I21|I22|I23|I241|I252", ukb_mi$f.41204.0.1), 1, 0)

ukb_mi$MI <- ukb_mi$primary + ukb_mi$primary_v2 + ukb_mi$secondary + ukb_mi$secondary_v2
ukb_mi$MI <- ifelse(ukb_mi$MI >= 1, 1, 0)

## Match with existing phenotype file
pheno <- left_join(ukb_prs, ukb_mi[, c("f.eid", "MI")], by = c("ID" = "f.eid"))

## Export
write.table(pheno, "ukbb_prs_pheno_w_MI.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

## Make GCTA pheno file
covar <- fread("build_covar_gcta.txt")
mi <- ukbb[ukbb$ID %in% covar$V1, c("ID", "ID", "MI")]
write.table(mi, "MI_v2/mi_build_pheno_gcta.txt", sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)

## Make PRS pheno files
build <- ukbb[ukbb$ID %in% covar$V1, c(1, 1:14, 32)]
test <- ukbb[!(ukbb$ID %in% covar$V1), c(1, 1:14, 32)]

write.table(build, "PRS/PLINK_MI/mi_pheno_build.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(test, "PRS/PLINK_MI/mi_pheno_test.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

############################################################
#################### Run GCTA Analysis #####################
############################################################

#!/usr/bin/env bash                                                                                                                    

#SBATCH --array=1-22%10                                                                                                                
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000                                                                            
#SBATCH --job-name=mi_gcta                                                                                                             
#printenv                                                                                                                              
#exit                                                                                                                                  

module load pkgsrc/2020Q3

gcta64 \
        --pfile LD_PLINK/ldsub_chr${SLURM_ARRAY_TASK_ID} \
        --fastGWA-lr \
        --pheno inc_prs/MI_v2/mi_build_pheno_gcta.txt \
        --qcovar inc_prs/build_qcovar_gcta.txt \
        --covar inc_prs/build_covar_gcta.txt \
        --threads 10 \
        --out inc_prs/MI_v2/mi_chr${SLURM_ARRAY_TASK_ID}
        
############################################################
######### Step 1 - Format FastGWA results into .ma #########
############################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

# files <- list.files(pattern=".SAIGE.bgen.txt")
files <- list.files(path = "inc_prs/MI_v2/", pattern=".fastGWA")

full <- data.table()
for(i in files){
	temp <- fread(paste0("inc_prs/MI_v2/", i))
	full <- rbind(full, temp)
}

full <- full[, c("SNP", "A1", "A2", "AF1", "BETA", "SE", "P", "N")]

# -------------------------
# Exclude on variation in N
# -------------------------

# Have a look at histogram
# hist(full$N,    xlab = "N - Per SNP sample size", main = "Locke et al.")

olr.lw <- median(full$N) - 1.5 * IQR(full$N)
olr.up <- median(full$N) + 1.5 * IQR(full$N)

# abline(v = olr.lw, col = 'red')
# abline(v = olr.up, col = 'red')

full.qc <- full[which(full$N > olr.lw & full$N  < olr.up), ] # 287105 variants lost
sd(full.qc$N)

write.table(full.qc, "mi_qc_n.ma", col.names = T, row.names = F, sep = " ", quote = F)

for (i in seq(0.1, 0.8, 0.1))
{
  print(paste0("Doing ", i))
  full.qc.pval <- full.qc[which(full.qc$P < i), ] 
  write.table(full.qc.pval, paste0("inc_prs_gctb/gctb/MI/GCTA_v2/mi_qc_n_pval_",i, ".ma"),
              col.names = T, row.names = F, sep = " ", quote = F)
}

########################################
######### Step 2 - Run sBayesR #########
########################################

#!/bin/bash                                                                                                                             

#SBATCH --job-name=ldm_mi                                                                                                             
#SBATCH --mem-per-cpu=60000                                                                                                             
#source gctb/gctb_2.02_Linux gctb                                                                                       

mkdir MI_qc_n/
mv mi_qc_n.ma MI_qc_n/

gctb/gctb_2.02_Linux/gctb \
    --sbayes R \
    --mldm inc_prs_gctb/test_data/ukbEURu_hm3_sparse_mldm_list.txt \
    --gwas-summary inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n/mi_qc_n.ma \
    --out inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n/mi_qc_n \
    --chain-length 2000 --burn-in 1000 --out-freq 10

#################################

#!/usr/bin/env bash

#SBATCH --array=1-8%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=ldm_qc_mi

#printenv
#exit                                                                                                          

mkdir MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}
mv mi_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}.ma MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}

#source gctb/gctb_2.02_Linux gctb

gctb/gctb_2.02_Linux/gctb \
    --sbayes R \
    --mldm inc_prs_gctb/test_data/ukbEURu_hm3_sparse_mldm_list.txt \
    --gwas-summary inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/mi_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}.ma \
    --out inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/mi_qc_n_pval_0.${SLURM_ARRAY_TASK_ID} \
    --chain-length 2000 --burn-in 1000 --out-freq 10
    

###############################################################
######### Step 3 - Make params file for plink scoring #########
###############################################################

library(data.table)
library(dplyr)

x <- read.table("MI_qc_n/mi_qc_n.snpRes", header = T)
x.2 <- x[,c("Name", "A1", "A1Effect")]
write.table(x.2, "MI_qc_n/mi_qc_n.params", row.names = F, sep = " ", quote = F)

for(i in 1:8){
	x <- read.table(paste0("MI_qc_n_pval_0.", i, "/mi_qc_n_pval_0.", i, ".snpRes"), header = T)
	x.2 <- x[,c("Name", "A1", "A1Effect")]
	write.table(x.2, paste0("MI_qc_n_pval_0.", i, "/mi_qc_n_pval_0.", i, ".params"), row.names = F, sep = " ", quote = F)
}

	
############################################
######### Step 4 - Run plink score #########
############################################

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000
#SBATCH --job-name=score_mi

#printenv
#exit

mkdir MI_qc_n/score

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --bfile LD_PLINK/LD_BED/ldsub_chr${SLURM_ARRAY_TASK_ID} \
    --exclude LD_PLINK/LD_BED/bim${SLURM_ARRAY_TASK_ID}.dups \
    --score inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n/mi_qc_n.params \
    --out inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n/score/mi_chr${SLURM_ARRAY_TASK_ID}

#################################
#################################

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000
#SBATCH --job-name=qc_n_p2

#printenv
#exit

mkdir MI_qc_n_pval_0.2/score

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --bfile LD_PLINK/LD_BED/ldsub_chr${SLURM_ARRAY_TASK_ID} \
    --exclude LD_PLINK/LD_BED/bim${SLURM_ARRAY_TASK_ID}.dups \
    --score inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n_pval_0.2/mi_qc_n_pval_0.2.params \
    --out inc_prs_gctb/gctb/MI/GCTA_v2/MI_qc_n_pval_0.2/score/mi_chr${SLURM_ARRAY_TASK_ID}


###############################################
######### Step 5 - Calculate Pred R^2 #########
################# NUMERIC HDL #################
###############################################

library(data.table)
library(dplyr)

ukb <- read.table("inc_prs/PRS/PLINK_MI/mi_pheno_test.txt", header = FALSE)

resLocke <- data.frame(matrix(0, nrow = 9, ncol = 2))
qc <- c("", "_pval_0.1", "_pval_0.2", "_pval_0.3", "_pval_0.4", "_pval_0.5", "_pval_0.6", "_pval_0.7", "_pval_0.8")

for(k in 1:length(qc)){
	lky.chr1 <- read.table(paste0("MI_qc_n", qc[k], "/score/mi_chr1.profile"), header = T)
	lky.all  <- lky.chr1
	for (i in seq(2, 22)){
# 	for(i in c(2:5, 7:22)){
		print(paste("Doing", i))
		lky.chri <- read.table(paste0("MI_qc_n", qc[k], "/score/mi_chr", i, ".profile"), header = T)
		lky.all$SCORE <- lky.all$SCORE + lky.chri$SCORE
  	}
  	lky.ukb <- inner_join(ukb, lky.all, by = c("V1" = "FID"))
  	resLocke[k, 1] <- paste0("MI_qc_n", qc[k])
	resLocke[k, 2] <- summary(lm(lky.ukb$V16 ~ lky.ukb$SCORE))$r.squared
	print(resLocke)
}
colnames(resLocke) <- c("QC", "Pred_R2")

write.table(resLocke, "mi_v2_results.txt", col.names = T, row.names = F, quote = F, sep = " ")

