##############################################################################
## Title: UKBB MI phenotype using date of diagnosis & SAIGE analysis
## Author: Regina Manansala
## Date Created: 14-December-2020
## Date Modified: 15-March-2021
##############################################################################

library(data.table)
library(dplyr)

## Import phenotype files
ukb_prs <- fread("inc_prs/ukbb_prs_pheno_w_bccount.txt")
test_ids <- fread("inc_prs_gctb/gctb/test_subset_IDs.txt") 
build <- fread("inc_prs/build_covar_gcta.txt")

## Import UKBB data
source('/ukb23544_v2.r')
inc <- bd
# source('/ukb25323_v2.r')
# gen <- bd
# source('/ukb29565_v2.r')
# edu <- bd
# source('/ukb40060.r')
# smk <- bd
# source("/ukb43095.r")
# dep <- bd
# source('/ukb21755_v2.r')

## Make MI variable
ukb_mi <- left_join(ukb_prs, inc[, c("f.eid", "f.42000.0.0")], by = c("ID"="f.eid"))
ukb_mi$mi_flag <- ifelse(is.na(ukb_mi$f.42000.0.0), 0, 1)

ukb_test <- subset(ukb_mi, ID %in% test_ids$V1)
ukb_build <- subset(ukb_mi, ID %in% build$V1)

table(ukb_test$mi_flag)
table(ukb_build$mi_flag)

############################################################
####################### SAIGE STEP 1 #######################
############################################################

#!/usr/bin/env bash

#SBATCH --array=15-16%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000
#SBATCH --job-name=step1

#printenv
#exit

source ~/Data/conda_install/bin/activate RSAIGE
    
Rscript ~/Data/manansa2/stroke/SAIGE/step1_fitNULLGLMM.R \
	--plinkFile=LD_PLINK/ldsub_chr${SLURM_ARRAY_TASK_ID} \
	--phenoFile=inc_prs_gctb/gctb/MI/SAIGE/ukbb_prs_pheno_mi_build.txt \
	--phenoCol=mi_flag \
	--covarColList=Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	--sampleIDColinphenoFile=ID \
	--traitType=binary \
	--outputPrefix=inc_prs_gctb/gctb/MI/SAIGE/mi \
	--nThreads=${SLURM_CPUS_PER_TASK} \
	--LOCO=FALSE \
	--IsOverwriteVarianceRatioFile=FALSE ## v0.38. Whether to overwrite the variance ratio file if the file already exists

conda deactivate

############################################################
####################### SAIGE STEP 2 #######################
############################################################


#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=step2_v1

#printenv
#exit

source ~/Data/conda_install/bin/activate RSAIGE

Rscript ~/Data/manansa2/stroke/SAIGE/step2_SPAtests.R \
	--bgenFile=chr${SLURM_ARRAY_TASK_ID}_MAF001_MI.bgen \
        --minMAF=0.001 \
        --minMAC=15 \
        --sampleFile=stroke/MI_sampleList_for_SAIGE.txt \
        --GMMATmodelFile=inc_prs_gctb/gctb/MI/SAIGE/mi.rda \
        --varianceRatioFile=inc_prs_gctb/gctb/MI/SAIGE/mi.varianceRatio.txt \
        --SAIGEOutputFile=inc_prs_gctb/gctb/MI/SAIGE/mi_chr${SLURM_ARRAY_TASK_ID}.SAIGE.bgen.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE

conda deactivate

#################################