##############################################################################
## Title: UKBB T2D sBayesR workflow
## Author: Regina Manansala
## Date Created: 23-November-2020
## Date Modified: 03-December-2020
##############################################################################

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)
# library(tidyverse)

## Read phenotype file
ukb_prs <- fread("ukbb_prs_pheno_w_bccount.txt")

## Read UKBB files
source('/raid-05/SPH/pauer/UKBB_PH/ukb23544_v2.r')
inc <- bd
source("/raid-05/SPH/pauer/UKBB_PH/ukb43095.r")
dep <- bd
source('/raid-05/SPH/pauer/UKBB_PH/ukb21755_v2.r')

## Get diabetes variables
inc_sub <- subset(inc, f.eid %in% ukb_prs$ID) %>% select(., c("f.eid", "f.2443.0.0", "f.20002.0.0"))
dep_sub <- subset(dep, f.eid %in% ukb_prs$ID) %>% select(., c("f.eid", "f.41202.0.0"))

t2d <- left_join(inc_sub, dep_sub, by = "f.eid")

## Create diabetes phenotype
t2d$diabetes <- ifelse(t2d$f.2443.0.0 == "Yes", 1, 0)
t2d$diabetes_self_reported <- ifelse(t2d$f.20002.0.0 == 1223, 1, 0)
t2d$diabetes_ic10 <- ifelse(grepl("E11", t2d$f.41202.0.0), 1, 0)

ukb_prs_t2d <- left_join(ukb_prs, t2d, by = c("ID" = "f.eid"))

## Export
write.table(ukb_prs_t2d, "ukbb_prs_pheno_w_t2d.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

##################################################################

## Make GCTA phenotype file

test_ids <- fread("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/test_subset_IDs.txt") 
build <- fread("build_covar_gcta.txt")

### UPDATE DIABETES VARIABLE
*** t2d_build_pheno_gcta <- fread("ukbb_prs_pheno_w_t2d.txt") %>% subset(., ID %in% build$V1) %>% select(., c("ID", "ID", "diabetes"))

write.table(t2d_build_pheno_gcta, "T2D/t2d_build_pheno_gcta.txt", sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)

##################################################################

## Run gwas using GCTA

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000
#SBATCH --job-name=t2d_gcta

#printenv
#exit

module load /raid-01/UITS/bacon/Pkgsrc/pkg-2019Q4/etc/modulefiles/pkgsrc/2019Q4

gcta64 \
	--pfile /raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/ldsub_chr${SLURM_ARRAY_TASK_ID} \
	--fastGWA-lr \
	--pheno /raid-04/SPH/pauer/manansa2/inc_prs/T2D/t2d_build_pheno_gcta.txt \
	--qcovar /raid-04/SPH/pauer/manansa2/inc_prs/build_qcovar_gcta.txt \
	--covar /raid-04/SPH/pauer/manansa2/inc_prs/build_covar_gcta.txt \
	--threads 10 \
	--out /raid-04/SPH/pauer/manansa2/inc_prs/T2D/t2d_chr${SLURM_ARRAY_TASK_ID}
	
