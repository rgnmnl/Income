##############################################################################
## Title: UKBB Height Analysis
## Author: Regina Manansala
## Date Created: 15-February-2019
## Date Modified: 17-February-2019
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(MASS)

#Import UKBB Income data
source('../ukb23544.r')
inc <- bd
source('../ukb25323.r')

#Import ID exclusion list
# ids <- as.list(read.table("w40458_20181016.txt", header=FALSE))
# bd[bd$f.eid %in% ids$V1, 1:5]
# ukbb_inc <- subset(bd, !(f.eid %in% ids$V1))

#Merge UKBB data and INC data
# foo <- ukbb_inc
ukbb_inc <- left_join(bd, inc, by='f.eid')

#Rename Variables (Income Data)
names(ukbb_inc)[names(ukbb_inc) == "f.eid"] <- "ID"
names(ukbb_inc)[names(ukbb_inc) == "f.738.0.0"] <- "income"
names(ukbb_inc)[names(ukbb_inc) == "f.50.0.0"] <- "height" 
names(ukbb_inc)[names(ukbb_inc) == "f.31.0.0"] <- "sex" 
names(ukbb_inc)[names(ukbb_inc) == "f.22001.0.0"] <- "sex_genetic" 
names(ukbb_inc)[names(ukbb_inc) == "f.54.0.0"] <- "loc_code"
names(ukbb_inc)[names(ukbb_inc) == "f.34.0.0"] <- "birth_yr"
names(ukbb_inc)[names(ukbb_inc) == "f.52.0.0"] <- "birth_mo"
names(ukbb_inc)[names(ukbb_inc) == "f.6142.0.0"] <- "empl_status"
names(ukbb_inc)[names(ukbb_inc) == "f.21001.0.0"] <- "BMI"
names(ukbb_inc)[names(ukbb_inc) == "f.21003.0.0"] <- "Age"
names(ukbb_inc)[names(ukbb_inc) == "f.132.0.0"] <- "Job"
names(ukbb_inc)[names(ukbb_inc) == "f.22617.0.0"] <- "SOC_Code"
names(ukbb_inc)[14:23] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

ukbb_inc$diabetes <- ifelse(is.na(ukbb_inc$f.2443.0.0) |ukbb_inc$f.2443.0.0 == "No" | ukbb_inc$f.2443.0.0 == "Prefer not to answer" | ukbb_inc$f.2443.0.0 == "Do not know", 0, 1) #diabetes
ukbb_inc$asthma <- ifelse(ukbb_inc$f.3786.0.0 == -3 | ukbb_inc$f.3786.0.0 == -1 | is.na(ukbb_inc$f.3786.0.0), 0, 1)
ukbb_inc$mi_flag <- ifelse(is.na(ukbb_inc$f.42000.0.0), 0, 1) #myocardial infarction
ukbb_inc$cancer_dx_dt <- ifelse(is.na(ukbb_inc$f.40005.0.0), 0, 1) #cancer diagnosis date
ukbb_inc$cancer_sr <- ifelse(is.na(ukbb_inc$f.20001.0.0), 0, 1) #self-reported cancer

ukbb_inc$race[ukbb_inc$f.21000.0.0 == "British"] <- "British"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Irish"] <- "Irish"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Any other white background" | ukbb_inc$f.21000.0.0 == "White"] <- "Other White"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Chinese"] <- "Chinese"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Caribbean" | ukbb_inc$f.21000.0.0 == "African" | ukbb_inc$f.21000.0.0 == "Any other Black background" | ukbb_inc$f.21000.0.0 == "Black or Black British"] <- "Black"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Mixed" | ukbb_inc$f.21000.0.0 == "White and Black Caribbean" | ukbb_inc$f.21000.0.0 == "White and Black African" | ukbb_inc$f.21000.0.0 == "White and Asian" | ukbb_inc$f.21000.0.0 == "Any other mixed background" | ukbb_inc$f.21000.0.0 == "Other ethnic group" | ukbb_inc$f.21000.0.0 == "Do not know" | ukbb_inc$f.21000.0.0 == "Prefer not to answer" | ukbb_inc$f.21000.0.0 == "Any other Asian background" | ukbb_inc$f.21000.0.0 == "Asian or Asian British"] <- "Other"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Indian" | ukbb_inc$f.21000.0.0 == "Pakistani" | ukbb_inc$f.21000.0.0 == "Bangladeshi"] <- "South Asian"

ukbb_inc_sub <- dplyr::select(ukbb_inc, -starts_with("f."))

#Remove NA income
inc_nomiss_inc <- subset(ukbb_inc_sub, !is.na(ukbb_inc_sub$income) & ukbb_inc_sub$income != "Prefer not to answer" & ukbb_inc_sub$income != "Do not know")
inc_nomiss_inc$income <- factor(inc_nomiss_inc$income)

#Only keep non-missing employment status
inc_nomiss_emp <- subset(inc_nomiss_inc, !(inc_nomiss_inc$empl_status %in% c("None of the above", "Prefer not to answer")))
inc_nomiss_emp$empl_status <- factor(inc_nomiss_emp$empl_status)

#Remove locations with <1000
inc_locex <- subset(inc_nomiss_emp, loc_code != 11023)

#Rename locations
inc_locex$location[inc_locex$loc_code == 11012] <- "Barts"
inc_locex$location[inc_locex$loc_code == 11021] <- "Birmingham"
inc_locex$location[inc_locex$loc_code == 11011] <- "Bristol"
inc_locex$location[inc_locex$loc_code == 11008] <- "Bury"
inc_locex$location[inc_locex$loc_code == 11003] <- "Cardiff"
inc_locex$location[inc_locex$loc_code == 11024] <- "Cheadle"
inc_locex$location[inc_locex$loc_code == 11020] <- "Croydon"
inc_locex$location[inc_locex$loc_code == 11005] <- "Edinburgh"
inc_locex$location[inc_locex$loc_code == 11004] <- "Glasgow"
inc_locex$location[inc_locex$loc_code == 11018] <- "Hounslow"
inc_locex$location[inc_locex$loc_code == 11010] <- "Leeds"
inc_locex$location[inc_locex$loc_code == 11016] <- "Liverpool"
inc_locex$location[inc_locex$loc_code == 11001] <- "Manchester"
inc_locex$location[inc_locex$loc_code == 11017] <- "Middlesborough"
inc_locex$location[inc_locex$loc_code == 11009] <- "Newcastle"
inc_locex$location[inc_locex$loc_code == 11013] <- "Nottingham"
inc_locex$location[inc_locex$loc_code == 11002] <- "Oxford"
inc_locex$location[inc_locex$loc_code == 11007] <- "Reading"
inc_locex$location[inc_locex$loc_code == 11014] <- "Sheffield"
inc_locex$location[inc_locex$loc_code == 11006] <- "Stoke"
inc_locex$location[inc_locex$loc_code == 11022] <- "Swansea"
# inc_locex$location[inc_locex$loc_code == 11023] <- "Wrexham"
inc_locex$location[inc_locex$loc_code == 11025] <- "Cheadle (Imaging)"
inc_locex$location[inc_locex$loc_code == 11026] <- "Reading (Imaging)"
inc_locex$location[inc_locex$loc_code == 11027] <- "Newcastle (imaging)"


#Remove NA race
inc_raceex <- subset(inc_locex, !is.na(race))
inc_whiteonly <- subset(inc_locex, race == "British")


#Limit Age to 25-65
inc_agesub <- subset(inc_whiteonly, Age >= 25 & Age <= 65)

#Create inc_analysis dataframe - remove missing PCs
inc_analysis <- subset(inc_agesub, !is.na(PC1))

#Create continuous income variable
inc_analysis$inc_con[inc_analysis$income == 'Less than 18,000'] <- 1
inc_analysis$inc_con[inc_analysis$income == '18,000 to 30,999'] <- 2
inc_analysis$inc_con[inc_analysis$income == '31,000 to 51,999'] <- 3
inc_analysis$inc_con[inc_analysis$income == '52,000 to 100,000'] <- 4
inc_analysis$inc_con[inc_analysis$income == 'Greater than 100,000'] <- 5

########

inc_analysis$race_fct <- as.factor(inc_analysis$race)

inc_analysis$employment <- as.character(inc_analysis$empl_status)
inc_analysis$employment <- as.factor(inc_analysis$employment)

inc_analysis$sexfct <- as.character(inc_analysis$sex)
inc_analysis$sexfct <- as.factor(inc_analysis$sexfct)

inc_analysis$location <- as.factor(inc_analysis$location)

#Double-check variables to exclude NAs/unwanted categories

# table(inc_analysis$income, exclude=NULL)
table(inc_analysis$employment, exclude=NULL)
table(inc_analysis$sexfct, exclude=NULL)
table(inc_analysis$race_fct, exclude=NULL)
table(inc_analysis$location, exclude=NULL)


########## CREATE COVARIATE SUBSET FOR HEIGHT ANALYSIS #################

#Subset columns and exclude NA Height
height_cov <- inc_analysis[!is.na(inc_analysis$height),c(1,1,2,22,34,30,4:13)]
names(height_cov)[1:6] <- c("#FID", "IID", "Height", "Age", "Sex", "Location")

#Random Sample 80% of data for analysis/20% for replication
# install.packages("ISLR")
library(ISLR)
smp_siz = floor(0.8*nrow(height_cov))
set.seed(123)
sample <- sample(seq_len(nrow(height_cov)),size = smp_siz)
ht_analysis <- height_cov[sample,]
ht_replicate <- height_cov[-sample,]

#Create Phenotype Files
write.table(ht_analysis, "HEIGHT_ANALYSIS/PLINK/ht_analysis.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(ht_replicate, "HEIGHT_ANALYSIS/PLINK/ht_replicate.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

#Stratify by Gender
write.table(ht_analysis[ht_analysis$Sex == "Female", -5], "HEIGHT_ANALYSIS/PLINK/ht_female.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(ht_analysis[ht_analysis$Sex == "Male", -5], "HEIGHT_ANALYSIS/PLINK/ht_male.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

write.table(ht_replicate[ht_replicate$Sex == "Female", -5], "HEIGHT_ANALYSIS/PLINK/ht_female_rep.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(ht_replicate[ht_replicate$Sex == "Male", -5], "HEIGHT_ANALYSIS/PLINK/ht_male_rep.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)


#########################################################
#########################################################
#########################################################

####SHELL SCRIPT - FULL GWAS

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=plink2_ht

#printenv
#exit

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink2 \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --pfile ~/Data/UKBB_Genotypes/chr${SLURM_ARRAY_TASK_ID}_PLINK_v2 \
    --out ~/Data/UKBB_PH/INC/APOE_CAREGIVER/inc_chr${SLURM_ARRAY_TASK_ID}_assoc \
    --glm hide-covar cols=chrom,pos,ref,alt1,a1count,totallele,a1freq,machr2,nobs,orbeta,se,tz,p \
    --pheno ~/Data/UKBB_PH/INC/HEIGHT_ANALYSIS/PLINK/ht_analysis.txt \
    --pheno-name Height \
    --covar ~/Data/UKBB_PH/INC/HEIGHT_ANALYSIS/PLINK/ht_analysis.txt \
    --covar-name Age,Sex,Location,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

#### SHELL SCRIPT - STRAT

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=plink2_ht

#printenv
#exit

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink2 \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --pfile ~/Data/UKBB_Genotypes/chr${SLURM_ARRAY_TASK_ID}_PLINK_v2 \
    --out ~/Data/UKBB_PH/INC/HEIGHT_ANALYSIS/PLINK/STRAT_MALE/inc_chr${SLURM_ARRAY_TASK_ID}_assoc \
    --glm hide-covar cols=chrom,pos,ref,alt1,a1count,totallele,a1freq,machr2,nobs,orbeta,se,tz,p \
    --pheno ~/Data/UKBB_PH/INC/HEIGHT_ANALYSIS/PLINK/ht_male.txt \
    --pheno-name Height \
    --covar ~/Data/UKBB_PH/INC/HEIGHT_ANALYSIS/PLINK/ht_male.txt \
    --covar-name Age,Location,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10