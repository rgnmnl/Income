##############################################################################
## Title: UKBB MI sBayesR Workflow
## Version: 2
## Note: For SAIGE gwas results (used all 150k samples)
## Author: Regina Manansala
## Date Created: 22-January-2021
## Date Modified: 02-February-2021
##############################################################################

## Run SAIGE Step 1

#!/usr/bin/env bash

#SBATCH --array=6%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=step1

#printenv
#exit

source ~/Data/conda_install/bin/activate RSAIGE
    
Rscript SAIGE/step1_fitNULLGLMM.R \
	--plinkFile=ldsub_chr${SLURM_ARRAY_TASK_ID} \
	--phenoFile=gctb/MI/SAIGE/ukbb_prs_pheno_mi_build.txt \
	--phenoCol=mi_flag \
	--covarColList=Age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	--sampleIDColinphenoFile=ID \
	--traitType=binary \
	--outputPrefix=gctb/MI/SAIGE/mi_chr${SLURM_ARRAY_TASK_ID} \
	--nThreads=${SLURM_CPUS_PER_TASK} \
	--LOCO=FALSE \
	--IsOverwriteVarianceRatioFile=TRUE ## v0.38. Whether to overwrite the variance ratio file if the file already exists

conda deactivate

## Run SAIGE Step 2

#!/usr/bin/env bash                                                                                                                     

#SBATCH --array=1-22%10                                                                                                                 
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000                                                                             
#SBATCH --job-name=step2_v1                                                                                                             

#printenv                                                                                                                               
#exit                                                                                                                                   

source ~/Data/conda_install/bin/activate RSAIGE

Rscript SAIGE/step2_SPAtests.R \
        --bgenFile=chr${SLURM_ARRAY_TASK_ID}_MAF001_MI_build.bgen \
        --minMAF=0.001 \
        --minMAC=15 \
        --sampleFile=MI_build_sampleList_for_SAIGE.txt \
        --GMMATmodelFile=gctb/MI/SAIGE/mi_chr${SLURM_ARRAY_TASK_ID}.rda \
        --varianceRatioFile=/raid-04/SPH/pauer/manansa2/inc_prs/gctb/MI/SAIGE/mi_chr${SLURM_ARRAY_TASK_ID}.varianceRatio.txt \
        --SAIGEOutputFile=gctb/MI/SAIGE/mi_chr${SLURM_ARRAY_TASK_ID}.SAIGE.bgen.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE

conda deactivate


### SBAYESR

############################################################
######### Step 1 - Format FastGWA results into .ma #########
############################################################

## NOTE: Used SNP information from LDL fastGWA results because sBayesR would not work with SAIGE results

#### https://www.dropbox.com/sh/9aougeoxw4ygo8k/AAD6PT3a3ggv1-KYHytbeUNha?dl=0&preview=ht_locke_heigwbc_wood_process.R

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

files <- list.files(path="gctb/MI/SAIGE/", pattern=".SAIGE.bgen.txt")

full <- data.table()
for(i in files){
	temp <- fread(paste0("gctb/MI/SAIGE/", i))
	full <- rbind(full, temp)
}

ldl_files <- list.files(path="LDL", pattern= ".fastGWA")
ldl <- data.table()
for(i in ldl_files){
	temp <- fread(paste0("LDL/", i))
	ldl <- rbind(ldl, temp)
}

mi_ldl <- left_join(ldl, full[, c("CHR", "POS", "rsid", "Allele1", "Allele2", "BETA", "SE", "p.value")], by = c("SNP"="rsid", "A1"="Allele1", "A2"="Allele2"))

# full$N_v2 <- round(full$N * full$imputationInfo,0)
# full$MAF <- ifelse(full$AF_Allele2 < 0.5, full$AF_Allele2, 1-full$AF_Allele2)
# full_maf <- subset(full, MAF > .01)

full_maf <- mi_ldl[, c("SNP", "A1", "A2", "AF1", "BETA.y", "SE.y", "p.value", "N")]

# -------------------------
# Exclude on variation in N
# -------------------------

# Have a look at histogram
# hist(full$N,    xlab = "N - Per SNP sample size", main = "Locke et al.")

olr.lw <- median(full_maf$N) - 1.5 * IQR(full_maf$N)
olr.up <- median(full_maf$N) + 1.5 * IQR(full_maf$N)

full.qc <- full_maf[which(full_maf$N > olr.lw & full_maf$N < olr.up), ] # 287105 variants lost

names(full.qc) <- c("SNP", "A1", "A2", "AF1", "BETA", "SE", "P", "N")

write.table(full.qc, "gctb/MI/MI_qc_n/mi_qc_n.ma", col.names = T, row.names = F, sep = " ", quote = F)

for (i in seq(0.1, 0.8, 0.1))
{
  print(paste0("Doing ", i))
  full.qc.pval <- full.qc[which(full.qc$P < i), ] 
#   print(dim(full.qc.pval))
  write.table(full.qc.pval, paste0("gctb/MI/MI_qc_n_pval_", i, "/mi_qc_n_pval_",i, ".ma"),
              col.names = T, row.names = F, sep = " ", quote = F)
}

###########################################################

# qc <- c("_pval_0.1", "_pval_0.2", "_pval_0.3", "_pval_0.4", "_pval_0.5", "_pval_0.6", "_pval_0.7", "_pval_0.8")
# 
# for(i in qc){
# 	temp <- fread(paste0("../MI_qc_n", i, "/mi_qc_n", i, ".ma"))
# 	names(temp) <- c("SNP", "A1", "A2", "MAF", "N", "BETA", "SE", "P")
# 	write.table(temp, paste0("../MI_qc_n", i, "/mi_qc_n", i, ".ma"), col.names = T, row.names = F, sep = " ", quote = F)
# }

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
    --gwas-summary gctb/MI/MI_qc_n/mi_qc_n.ma \
    --out gctb/MI/MI_qc_n/mi_qc_n \
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
    --gwas-summary gctb/MI/MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/mi_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}.ma \
    --out gctb/MI/MI_qc_n_pval_0.${SLURM_ARRAY_TASK_ID}/mi_qc_n_pval_0.${SLURM_ARRAY_TASK_ID} \
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
#SBATCH --job-name=score_test

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
    --score gctb/MI/MI_qc_n/mi_qc_n.params \
    --out gctb/MI/MI_qc_n/score/mi_chr${SLURM_ARRAY_TASK_ID}
    
    

###############################################
######### Step 5 - Calculate Pred R^2 #########
################# NUMERIC HDL #################
###############################################

library(data.table)
library(dplyr)
library(rcompanion)

ukb <- read.table("gctb/MI/SAIGE/ukbb_prs_pheno_mi_test.txt", header = TRUE)

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
  	lky.ukb <- inner_join(ukb, lky.all, by = c("ID" = "FID"))
  	resLocke[k, 1] <- paste0("MI_qc_n", qc[k])
  	mod <- glm(lky.ukb$mi_flag ~ lky.ukb$SCORE, family="binomial")
	resLocke[k, 2] <- nagelkerke(mod)$Pseudo.R.squared.for.model.vs.null[3]
	print(resLocke)
}
colnames(resLocke) <- c("QC", "Pred_R2")

write.table(resLocke, "mi_saige_results_v2.txt", col.names = T, row.names = F, quote = F, sep = " ")
























