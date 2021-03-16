##############################################################################
## Title: PRS Genotype File Creation
## Author: Regina Manansala
## Date Created: 27-December-2019
## Date Modified: 30-October-2020
##############################################################################

#########################################
############ CREATE SNP LIST ############
#########################################

### Mortimer Path: inc_prs/SNP_List

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)
library(tidyverse)
#library(nnet)

## Read SNP files
inc_snps <- fread("Income_SNPs.txt")
ht_snps <- fread("Height_SNPs.txt")
ed_snps <- fread("EduYears_Lead_SNPs.txt")
smk_snps <- fread("Smoking_SNPs.txt")
alc_snps <- fread("Alcohol_SNPs.txt")
dep_snps <- fread("Depression_SNPs.txt")

## Combine all SNPs
inc_snps1 <- inc_snps %>% mutate(Phenotype = "Income", SNP_1 = paste(inc_snps$Independent_significant_SNPs, inc_snps$effect_allele, sep = "_"), SNP_2 = paste(inc_snps$Independent_significant_SNPs, inc_snps$non_effect_allele, sep = "_")) %>% rename(.,SNP = Independent_significant_SNPs) %>% select(., c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome"))
inc_snps2 <- inc_snps %>% mutate(Phenotype = "Income", SNP_1 = paste(inc_snps$Tagged_SNPs, inc_snps$effect_allele, sep = "_"), SNP_2 = paste(inc_snps$Tagged_SNPs, inc_snps$non_effect_allele, sep = "_")) %>% rename(.,SNP = Tagged_SNPs) %>% select(., c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome"))

inc_snps3 <- inc_snps1[inc_snps1$SNP_1 %in% inc_snps2$SNP_1 & inc_snps1$SNP_2 %in% inc_snps2$SNP_2,] %>% subset(., !duplicated(SNP_1, SNP_2)) 
ht_snps <- ht_snps %>% mutate(Phenotype = "Height", SNP_1 = paste(ht_snps$SNP, ht_snps$Tested_Allele, sep = "_"), SNP_2 = paste(ht_snps$SNP, ht_snps$Other_Allele, sep = "_")) %>% rename(., Chromosome = CHR) %>% select(., c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome"))
ed_snps <- ed_snps %>% mutate(Phenotype = "Education", SNP_1 = paste(ed_snps$SNP, ed_snps$`Allele 1`, sep = "_"), SNP_2 = paste(ed_snps$SNP, ed_snps$Allele2, sep = "_")) %>% rename(., Chromosome = Chr) %>% select(., c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome"))
smk_snps <- smk_snps %>% mutate(Phenotype = "Smoking", SNP_1 = paste(smk_snps$rsID, smk_snps$Reference_Allele, sep = "_"), SNP_2 = paste(smk_snps$rsID, smk_snps$Alternate_Allele, sep = "_")) %>% rename(.,SNP = rsID, Chromosome = Chr) %>% select(., c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome"))
alc_snps <- alc_snps %>% mutate(Phenotype = "Alcohol", SNP_1 = paste(alc_snps$rsID, alc_snps$Reference_Allele, sep = "_"), SNP_2 = paste(alc_snps$rsID, alc_snps$Alternate_Allele, sep = "_")) %>% rename(.,SNP = rsID, Chromosome = Chr) %>% select(., c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome"))
dep_snps <- dep_snps %>% mutate(Phenotype = "Depression", SNP_1 = paste(dep_snps$SNP, dep_snps$Tested_Allele, sep = "_"), SNP_2 = paste(dep_snps$SNP, dep_snps$A2, sep = "_")) %>% rename(., Chromosome = Chr) %>% select(., c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome"))

tyrrell <- fread("/raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES/Tyrrell_2016_BMI_SNPS.txt", header = F)
tyrrell_pleio <- fread("/raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES/Tyrrell_2016_BMI_Pleio_SNPS.txt", header = F)
speliotes <- fread("/raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES/Speliotes_2010_BMI_SNPs.txt", header = F)
locke <- fread("/raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES/Locke_2015_BMI_SNPs.txt", header = F)
locke_pleio <- fread("/raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES/Locke_Pleio_SNPs.txt", header = T)

union_bmi_imp <- fread("/raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES/union_bmi_imp.txt")
bmi_snps <- union_bmi_imp[,grep("rs", colnames(union_bmi_imp), value = TRUE)] %>% as.data.frame()
bmi_snps <- cbind("BMI", bmi_snps, NA, NA, NA)
names(bmi_snps) <- c("Phenotype", "SNP", "SNP_1", "SNP_2", "Chromosome")

all_snps <- rbind(inc_snps3, ht_snps, ed_snps, smk_snps, alc_snps, dep_snps, bmi_snps)

all_snps <- rbind(inc_snps3[, c("Phenotype", "SNP", "SNP_1", "SNP_2")], ht_snps[, c("Phenotype", "SNP", "SNP_1", "SNP_2")], ed_snps[, c("Phenotype", "SNP", "SNP_1", "SNP_2")], smk_snps[, c("Phenotype", "SNP", "SNP_1", "SNP_2")], alc_snps[, c("Phenotype", "SNP", "SNP_1", "SNP_2")], dep_snps[, c("Phenotype", "SNP", "SNP_1", "SNP_2")])
# dim(all_snps[duplicated(all_snps$SNP_1) | duplicated(all_snps$SNP_1, fromLast = TRUE) | duplicated(all_snps$SNP_2) | duplicated(all_snps$SNP_2, fromLast = TRUE), ])
# all_snps <- all_snps %>% arrange(SNP) %>% subset(., !duplicated(all_snps$SNP_1) | !duplicated(all_snps$SNP_2))
write.table(all_snps, "all_snps_w_chr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#############################################
############ CONVERT PGEN TO BED ############
#############################################

#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=48000
#SBATCH --job-name=pgen_to_bed

#printenv
#exit

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink2 \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --pfile ~/Data/UKBB_Genotypes/chr${SLURM_ARRAY_TASK_ID}_PLINK_v2 \
    #--extract ~/Data/UKBB_PH/INC/GENO_FILES_12262019/ht_ed_inc_snps.txt \
    --out ~/Data/UKBB_PH/v2_BED/chr${SLURM_ARRAY_TASK_ID}_v2 \
    --make-bed


############################################
############ LD Pruning OG Data ############
############################################
#!/usr/bin/env bash

#SBATCH --array=1-22%10
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=60000
#SBATCH --job-name=ld_subset

#printenv
#exit

source /etc/bashrc
module purge
module load pkgsrc/2018Q4

plink2 \
    --memory ${SLURM_MEM_PER_NODE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --bfile ~/raid-05/SPH/pauer/UKBB/UKBB_G_Income/chr${SLURM_ARRAY_TASK_ID}_v2 \
    --indep-pairwise 500 50 0.1 \
    --out ~/Data/manansa2/inc_prs/PLINK/ldsub_chr${SLURM_ARRAY_TASK_ID} \
    --make-bed

##############################################
############## CREATE MR FILES ###############
##############################################

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)
library(tibble)

all_snps <- fread("all_snps.txt")

geno_path <- "/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK"
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

wg <- ColumnLinkedMatrix(chr1, chr2, chr3, chr4, chr5, chr6, chr7,
						chr8, chr9, chr10, chr11, chr12, chr13, chr14,
						chr15, chr16, chr17, chr18, chr19, chr20, chr21)
geno <- wg[, colnames(wg) %in% c(all_snps$SNP_1, all_snps$SNP_2)]
#write.table(geno, "geno.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
wg_sub <- geno %>% as.data.frame() #%>% rownames_to_column('ID') %>% mutate(ID = gsub("_.*","",ID))
wg_sub$ID <- gsub("_.*","",rownames(wg_sub)) %>% as.integer()

wg_sub_imp <- do.call(cbind, lapply(wg_sub[,2:8528], function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x}))

dim(wg[, colnames(wg) %in% all_snps$V1])						
geno <- as.matrix(wg) %>% as.data.frame()
geno$ID <- gsub("_.*","",rownames(geno)) %>% as.integer()
foo <- gsub("_.*", "", colnames(geno))
foo[43] <- paste0(foo[43], ".1")
names(geno) <- foo

wg_sub <- fread("ld_geno.txt")

############## Get Phenotypes ###############

source('/raid-05/SPH/pauer/UKBB_PH/ukb23544_v2.r')
inc <- bd
source('/raid-05/SPH/pauer/UKBB_PH/ukb25323_v2.r')
gen <- bd
source('/raid-05/SPH/pauer/UKBB_PH/ukb29565_v2.r')
edu <- bd
source('/raid-05/SPH/pauer/UKBB_PH/ukb40060.r')
smk <- bd

ukbb <- left_join(smk, edu, by='f.eid') %>% left_join(gen, by = 'f.eid') %>% left_join(inc, by = 'f.eid')

#write.table(ukbb, "ukbb.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

ukbb$income[ukbb$f.738.0.0 == 'Less than 18,000'] <- 1
ukbb$income[ukbb$f.738.0.0 == '18,000 to 30,999'] <- 2
ukbb$income[ukbb$f.738.0.0 == '31,000 to 51,999'] <- 3
ukbb$income[ukbb$f.738.0.0 == '52,000 to 100,000'] <- 4
ukbb$income[ukbb$f.738.0.0 == 'Greater than 100,000'] <- 5

names(ukbb)[names(ukbb) == "f.eid"] <- "ID"
#names(ukbb)[names(ukbb) == "f.738.0.0"] <- "income"
names(ukbb)[names(ukbb) == "f.50.0.0"] <- "height" 
names(ukbb)[names(ukbb) == "f.31.0.0"] <- "sex" 
names(ukbb)[names(ukbb) == "f.22001.0.0"] <- "sex_genetic" 
names(ukbb)[names(ukbb) == "f.54.0.0"] <- "loc_code"
names(ukbb)[names(ukbb) == "f.34.0.0"] <- "birth_yr"
names(ukbb)[names(ukbb) == "f.52.0.0"] <- "birth_mo"
names(ukbb)[names(ukbb) == "f.6142.0.0"] <- "empl_status"
names(ukbb)[names(ukbb) == "f.21001.0.0"] <- "BMI"
names(ukbb)[names(ukbb) == "f.21003.0.0"] <- "Age"
names(ukbb)[names(ukbb) == "f.132.0.0"] <- "Job"
names(ukbb)[names(ukbb) == "f.22617.0.0"] <- "SOC_Code"
names(ukbb)[names(ukbb) == "f.1558.0.0"] <- "alc_int_freq"
names(ukbb)[names(ukbb) == "f.20160.0.0"] <- "ever_smoked"
names(ukbb)[names(ukbb) == "f.20414.0.0"] <- "alc_dr_freq"
names(ukbb)[names(ukbb) == "f.20544.0.1"] <- "mental_health"
names(ukbb)[names(ukbb) == "f.845.0.0"] <- "Age_Ed"

names(ukbb)[42:51] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

ukbb$diabetes <- ifelse(is.na(ukbb$f.2443.0.0) |ukbb$f.2443.0.0 == "No" | ukbb$f.2443.0.0 == "Prefer not to answer" | ukbb$f.2443.0.0 == "Do not know", 0, 1) #diabetes
ukbb$asthma <- ifelse(ukbb$f.3786.0.0 == -3 | ukbb$f.3786.0.0 == -1 | is.na(ukbb$f.3786.0.0), 0, 1)
ukbb$mi_flag <- ifelse(is.na(ukbb$f.42000.0.0), 0, 1) #myocardial infarction
ukbb$cancer_dx_dt <- ifelse(is.na(ukbb$f.40005.0.0), 0, 1) #cancer diagnosis date
ukbb$cancer_sr <- ifelse(is.na(ukbb$f.20001.0.0), 0, 1) #self-reported cancer

ukbb$race[ukbb$f.21000.0.0 == "British"] <- "British"
ukbb$race[ukbb$f.21000.0.0 == "Irish"] <- "Irish"
ukbb$race[ukbb$f.21000.0.0 == "Any other white background" | ukbb$f.21000.0.0 == "White"] <- "Other White"
ukbb$race[ukbb$f.21000.0.0 == "Chinese"] <- "Chinese"
ukbb$race[ukbb$f.21000.0.0 == "Caribbean" | ukbb$f.21000.0.0 == "African" | ukbb$f.21000.0.0 == "Any other Black background" | ukbb$f.21000.0.0 == "Black or Black British"] <- "Black"
ukbb$race[ukbb$f.21000.0.0 == "Mixed" | ukbb$f.21000.0.0 == "White and Black Caribbean" | ukbb$f.21000.0.0 == "White and Black African" | ukbb$f.21000.0.0 == "White and Asian" | ukbb$f.21000.0.0 == "Any other mixed background" | ukbb$f.21000.0.0 == "Other ethnic group" | ukbb$f.21000.0.0 == "Do not know" | ukbb$f.21000.0.0 == "Prefer not to answer" | ukbb$f.21000.0.0 == "Any other Asian background" | ukbb$f.21000.0.0 == "Asian or Asian British"] <- "Other"
ukbb$race[ukbb$f.21000.0.0 == "Indian" | ukbb$f.21000.0.0 == "Pakistani" | ukbb$f.21000.0.0 == "Bangladeshi"] <- "South Asian"

ukbb$location[ukbb$loc_code == 11012] <- "Barts"
ukbb$location[ukbb$loc_code == 11021] <- "Birmingham"
ukbb$location[ukbb$loc_code == 11011] <- "Bristol"
ukbb$location[ukbb$loc_code == 11008] <- "Bury"
ukbb$location[ukbb$loc_code == 11003] <- "Cardiff"
ukbb$location[ukbb$loc_code == 11024] <- "Cheadle"
ukbb$location[ukbb$loc_code == 11020] <- "Croydon"
ukbb$location[ukbb$loc_code == 11005] <- "Edinburgh"
ukbb$location[ukbb$loc_code == 11004] <- "Glasgow"
ukbb$location[ukbb$loc_code == 11018] <- "Hounslow"
ukbb$location[ukbb$loc_code == 11010] <- "Leeds"
ukbb$location[ukbb$loc_code == 11016] <- "Liverpool"
ukbb$location[ukbb$loc_code == 11001] <- "Manchester"
ukbb$location[ukbb$loc_code == 11017] <- "Middlesborough"
ukbb$location[ukbb$loc_code == 11009] <- "Newcastle"
ukbb$location[ukbb$loc_code == 11013] <- "Nottingham"
ukbb$location[ukbb$loc_code == 11002] <- "Oxford"
ukbb$location[ukbb$loc_code == 11007] <- "Reading"
ukbb$location[ukbb$loc_code == 11014] <- "Sheffield"
ukbb$location[ukbb$loc_code == 11006] <- "Stoke"
ukbb$location[ukbb$loc_code == 11022] <- "Swansea"
# ukbb$location[ukbb$loc_code == 11023] <- "Wrexham"
ukbb$location[ukbb$loc_code == 11025] <- "Cheadle (Imaging)"
ukbb$location[ukbb$loc_code == 11026] <- "Reading (Imaging)"
ukbb$location[ukbb$loc_code == 11027] <- "Newcastle (imaging)"

ukbb_sub <- dplyr::select(ukbb, -starts_with("f."))

write.table(ukbb_sub, "ukbb_var_sub.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

tyrrel_bmi <- fread("/raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES/tyrrell_bmi.txt")
ukbb_id_sub <- ukbb_sub[ukbb_sub$ID %in% tyrrel_bmi$ID,]
ukbb_geno <- left_join(ukbb_id_sub, wg_sub, by = "ID")
ukbb_geno <- ukbb_geno[, c(1, 26, 19, 35, 9:18, 2, 3, 5, 7, 25, 28, 6, 36:8562)] %>% left_join(., tyrrel_bmi[, c(1, 18:93)], by = "ID")

foo <- ukbb_geno[rowSums(is.na(ukbb_geno[, 22:8548])) == 8527,]

geno_imp <- do.call(cbind, lapply(geno[,2:8528], function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x}))
geno_imp <- do.call(cbind, lapply(geno_sub, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x}))

imp <- as.data.frame(ukbb_nomiss)
for(i in 22:8548){
	imp[is.na(imp[, i]), i] <- mean(imp[,i], na.rm = TRUE)
}

############## Prelim Models ##############

imp <- fread("ukbb_prs.txt") %>% as.data.frame()
imp$alc[imp$alc_int_freq == "Daily or almost daily"] <- 7
imp$alc[imp$alc_int_freq == "Three or four times a week"] <- 4
imp$alc[imp$alc_int_freq == "Once or twice a week"] <- 2
imp$alc[imp$alc_int_freq == "One to three times a month"] <- .03
imp$alc[imp$alc_int_freq %in% c("Special occasions only", "Never", "Prefer not to answer")] <- 0
imp <- imp %>% mutate(smk = ifelse(ever_smoked == "Yes", 1, 0), depression = ifelse(mental_health != "Depression" | is.na(mental_health), 0, 1), alc = as.numeric(alc))

names(imp[!grepl("rs", names(imp))])

snp_list <- str_c(names(imp)[22:8548],collapse='+')
short_list <- str_c(names(imp)[22:30],collapse='+')

for(i in c("alc", "smk", "depression", "height", "income", "Age_Ed")){
print(length(table(imp[[i]])))
}

mod <- glm(ever_smoked ~ imp[,22:30], data=imp, family = binomial)

mod <- glm(paste0("smk ~ ",short_list), data=imp, family = binomial)

for(i in c("alc", "smk", "depression", "height", "income", "Age_Ed")){
	if(length(table(imp[[i]])) > 2){
		mod <- lm(paste0(i, " ~ ",short_list), data=imp)
		assign(i, summary(mod)$coef)
	} else {
		mod <- glm(paste0(i, " ~ ",short_list), data=imp, family = binomial)
		assign(i, summary(mod)$coef)
	}
}

for(i in 1:10){#nrow(M_val)
	x <- M_val[i,]
	mod <- lm(x ~ relevel(ph2$quin07di, ref = "0") + ph2$AgeMom + relevel(ph2$Eth9, ref = "1") + relevel(ph2$eduma, ref = "4") + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[[i]] <- summary(mod)$coef[2,]
	if(i %in% c(1, nrow(M_val))) {
		print(Sys.time())
		print("Model 6")
	}
}

names(results) <- rownames(M_val)
big_data <- do.call(rbind, results) %>% as.data.frame()
names(big_data)[1:4] <- c("beta", "se", "t", "pvalue")
