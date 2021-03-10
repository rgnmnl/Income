##############################################################################
## Title: UKBB Education sBayesR Workflow
## Version: 1
## Author: Regina Manansala
## Date Created: 15-September-2020
## Date Modified: 15-September-2020
##############################################################################

############################################################
######### Step 1 - Format FastGWA results into .ma #########
############################################################

#### https://www.dropbox.com/sh/9aougeoxw4ygo8k/AAD6PT3a3ggv1-KYHytbeUNha?dl=0&preview=bmi_locke_height_wood_process.R

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

files <- list.files(path="EDU", pattern= ".fastGWA")

full <- data.table()
for(i in files){
	temp <- fread(paste0("EDU/", i))
	full <- rbind(full, temp)
# 	assign(inc_gwa, full)
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

full.qc <- full[which(full$N > olr.lw & full$N  < olr.up), ] # 452717 variants lost
sd(full.qc$N)

write.table(full.qc, "edu_qc_n.ma", col.names = T, row.names = F, sep = " ", quote = F)

for (i in seq(0.1, 0.8, 0.1))
{
  print(paste0("Doing ", i))
  full.qc.pval <- full.qc[which(full.qc$P < i), ] 
  write.table(full.qc.pval, paste0("gctb/EDU/edu_qc_n_pval_",i, ".ma"),
              col.names = T, row.names = F, sep = " ", quote = F)
}

########################################
######### Step 2 - Run sBayesR #########
########################################

#!/bin/bash                                                                                                                             

#SBATCH --job-name=ldm_edu                                                                                                             
#SBATCH --mem-per-cpu=60000                                                                                                             
#source ~/Data/manansa2/gctb/gctb_2.02_Linux gctb                                                                                       

mkdir EDU_qc_n/

~/Data/manansa2/gctb/gctb_2.02_Linux/gctb \
    --sbayes R \
    --mldm test_data/ukbEURu_hm3_sparse_mldm_list.txt \
    --gwas-summary gctb/EDU/edu_qc_n.ma \
    --out gctb/EDU/EDU_qc_n/edu_qc_n \
    --chain-length 2000 --burn-in 1000 --out-freq 10
    
#################################

#!/usr/bin/env bash

#SBATCH --array=1-8%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=sbR_edu

#printenv
#exit                                                                                                          

#mkdir EDU_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}
mv edu_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}.ma EDU_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}

#source ~/Data/manansa2/gctb/gctb_2.02_Linux gctb

~/Data/manansa2/gctb/gctb_2.02_Linux/gctb \
    --sbayes R \
    --mldm test_data/ukbEURu_hm3_sparse_mldm_list.txt \
    --gwas-summary gctb/EDU/EDU_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/edu_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}.ma \
    --out gctb/EDU/EDU_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/edu_qc_n_pval_0.${SLURM_ARRAY_TASK_ID} \
    --chain-length 2000 --burn-in 1000 --out-freq 10
    
###############################################################
######### Step 3 - Make params file for plink scoring #########
###############################################################

library(data.table)
library(dplyr)

x <- read.table("EDU_qc_n/edu_qc_n.snpRes", header = T)
x.2 <- x[,c("Name", "A1", "A1Effect")]
write.table(x.2, "EDU_qc_n/edu_qc_n.params", row.names = F, sep = " ", quote = F)

for(i in 1:8){
	x <- read.table(paste0("EDU_qc_n_pval_0.", i, "/edu_qc_n_pval_0.", i, ".snpRes"), header = T)
	x.2 <- x[,c("Name", "A1", "A1Effect")]
	write.table(x.2, paste0("EDU_qc_n_pval_0.", i, "/edu_qc_n_pval_0.", i, ".params"), row.names = F, sep = " ", quote = F)
}

