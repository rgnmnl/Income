##############################################################################
## Title: UKBB MI sBayesR Workflow
## Author: Regina Manansala
## Date Created: 22-January-2021
## Date Modified: 02-February-2021
##############################################################################

############################################################
######### Step 1 - Format FastGWA results into .ma #########
############################################################

#### https://www.dropbox.com/sh/9aougeoxw4ygo8k/AAD6PT3a3ggv1-KYHytbeUNha?dl=0&preview=ht_locke_heigwbc_wood_process.R

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

# files <- list.files(pattern=".SAIGE.bgen.txt")
files <- list.files(path = "MI/", pattern=".fastGWA")

full <- data.table()
for(i in files){
	temp <- fread(paste0("MI/", i))
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
  write.table(full.qc.pval, paste0("gctb/MI/GCTA/mi_qc_n_pval_",i, ".ma"),
              col.names = T, row.names = F, sep = " ", quote = F)
}

########################################
######### Step 2 - Run sBayesR #########
########################################

#!/bin/bash                                                                                                                             

#SBATCH --job-name=ldm_mi                                                                                                             
#SBATCH --mem-per-cpu=60000                                                                                                             
#source ~/Data/manansa2/gctb/gctb_2.02_Linux gctb                                                                                       

mkdir MI_qc_n/
mv mi_qc_n.ma MI_qc_n/

~/Data/manansa2/gctb/gctb_2.02_Linux/gctb \
    --sbayes R \
    --mldm test_data/ukbEURu_hm3_sparse_mldm_list.txt \
    --gwas-summary gctb/MI/GCTA/MI_qc_n/mi_qc_n.ma \
    --out gctb/MI/GCTA/MI_qc_n/mi_qc_n \
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

#source ~/Data/manansa2/gctb/gctb_2.02_Linux gctb

~/Data/manansa2/gctb/gctb_2.02_Linux/gctb \
    --sbayes R \
    --mldm test_data/ukbEURu_hm3_sparse_mldm_list.txt \
    --gwas-summary gctb/MI/GCTA/MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/mi_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}.ma \
    --out gctb/MI/GCTA/MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/mi_qc_n_pval_0.${SLURM_ARRAY_TASK_ID} \
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
    --bfile ldsub_chr${SLURM_ARRAY_TASK_ID} \
    --exclude bim${SLURM_ARRAY_TASK_ID}.dups \
    --score gctb/MI/GCTA/MI_qc_n/mi_qc_n.params \
    --out gctb/MI/GCTA/MI_qc_n/score/mi_chr${SLURM_ARRAY_TASK_ID}

#################################
#################################

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000
#SBATCH --job-name=qc_n_p8

#printenv
#exit

mkdir MI_qc_n_pval_0.8/score

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --bfile ldsub_chr${SLURM_ARRAY_TASK_ID} \
    --exclude bim${SLURM_ARRAY_TASK_ID}.dups \
    --score gctb/MI/GCTA/MI_qc_n_pval_0.8/mi_qc_n_pval_0.8.params \
    --out gctb/MI/GCTA/MI_qc_n_pval_0.8/score/mi_chr${SLURM_ARRAY_TASK_ID}


###############################################
######### Step 5 - Calculate Pred R^2 #########
################# NUMERIC HDL #################
###############################################

library(data.table)
library(dplyr)

ukb <- read.table("MI/mi_test_20k_subset.txt", header = TRUE)

resLocke <- data.frame(matrix(0, nrow = 9, ncol = 2))
qc <- c("", "_pval_0.1", "_pval_0.2", "_pval_0.3", "_pval_0.4", "_pval_0.5", "_pval_0.6", "_pval_0.7", "_pval_0.8")

resLocke <- data.frame(matrix(0, nrow = 2, ncol = 2))
qc <- c("", "_pval_0.8")

for(k in 1:length(qc)){
	lky.chr1 <- read.table(paste0("MI_qc_n", qc[k], "/score/mi_chr1.profile"), header = T)
	lky.all  <- lky.chr1
	for (i in seq(2, 22)){
# 	for(i in c(2:5, 7:22)){
		print(paste("Doing", i))
		lky.chri <- read.table(paste0("MI_qc_n", qc[k], "/score/mi_chr", i, ".profile"), header = T)
		lky.all$SCORE <- lky.all$SCORE + lky.chri$SCORE
  	}
  	lky.ukb <- inner_join(ukb, lky.all, by = c("ID" = "FID"))
  	resLocke[k, 1] <- paste0("MI_qc_n", qc[k])
	resLocke[k, 2] <- summary(lm(lky.ukb$mi_flag ~ lky.ukb$SCORE))$r.squared
	print(resLocke)
}
colnames(resLocke) <- c("QC", "Pred_R2")

write.table(resLocke, "mi_20k_results.txt", col.names = T, row.names = F, quote = F, sep = " ")


	
	








