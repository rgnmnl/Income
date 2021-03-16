##############################################################################
## Title: UKBB MI GCTA
## Note: Run MI gwas using GCTA
## Author: Regina Manansala
## Date Created: 26-January-2021
## Date Modified: 26-January-2021
##############################################################################

library(data.table)
library(dplyr)

## Read phenotype files
mi <- fread("inc_prs_gctb/gctb/MI/SAIGE/ukbb_prs_pheno_mi.txt")
build <- fread("inc_prs_gctb/gctb/MI/SAIGE/ukbb_prs_pheno_mi_build.txt")
test <- fread("inc_prs_gctb/gctb/MI/SAIGE/ukbb_prs_pheno_mi_test.txt")

## Read original covariate files
covar <- fread("build_covar_gcta.txt")
qcovar <- fread("build_qcovar_gcta.txt")

## Create 20k subset
build_cases <- subset(build, mi_flag == 1)
build_controls <- subset(build, mi_flag == 0) %>% sample_n(., 15000)
build_mi <- bind_rows(build_cases, build_controls)

covar_sub <- subset(covar, V1 %in% build_mi$ID)
qcovar_sub <- subset(qcovar, V1 %in% build_mi$ID) 

test_cases <- subset(test, mi_flag == 1)
test_controls <- subset(test, mi_flag == 0) %>% sample_n(., 15000)
test_mi <- bind_rows(test_cases, test_controls)

## Export relevant GCTA files
write.table(covar_sub, "MI/mi_build_covar_gcta.txt", sep = "\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
write.table(qcovar_sub, "MI/mi_build_qcovar_gcta.txt", sep = "\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
write.table(trait, "MI/mi_build_pheno_gcta.txt", sep = "\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
write.table(build_mi, "MI/mi_build_20k_subset.txt", sep = "\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
write.table(test_mi, "MI/mi_test_20k_subset.txt", sep = "\t", col.names=TRUE, row.names=FALSE, quote = FALSE)

################################################################################

## Run GCTA

#!/usr/bin/env bash                                                                                                                     

#SBATCH --array=1-22%10                                                                                                                 
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000                                                                             
#SBATCH --job-name=mi_gcta                                                                                                             

#printenv                                                                                                                               
#exit                                                                                                                                   

module load /raid-01/UITS/bacon/Pkgsrc/pkg-2019Q4/etc/modulefiles/pkgsrc/2019Q4

gcta64 \
        --pfile ldsub_chr${SLURM_ARRAY_TASK_ID} \
        --fastGWA-lr \
        --pheno inc_prs/MI/mi_build_20k_subset.txt \
        --qcovar inc_prs/MI/mi_build_qcovar_gcta.txt \
        --covar inc_prs/MI/mi_build_covar_gcta.txt \
        --threads 10 \
        --out inc_prs/MI/mi_chr${SLURM_ARRAY_TASK_ID}