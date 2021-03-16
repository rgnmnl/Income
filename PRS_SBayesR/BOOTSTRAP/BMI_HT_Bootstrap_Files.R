##############################################################################
## Title: Make BMI, HT Bootstrap Files
## Author: Regina Manansala
## Date Created: 10-March-2021
## Date Modified: 10-March-2021
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)

ukbb <- fread("../../ukbb_prs_complete_v2.txt")
test_ids <- fread("inc_prs_gctb/gctb/test_subset_IDs.txt") 

# Create BMI subset
## use this when subsetting bootstrap samples within the loop
test <- ukbb[ukbb$ID %in% test_ids$V1, c(1, 18:19, 8549:8642)]
build <- ukbb[!(ukbb$ID %in% test_ids$V1), c(1, 18:19, 8549:8642)]

# Get list of BMI PRS files
bmi_files <- list.files(path="inc_prs_gctb/gctb/BMI/BMI_qc_n/score", pattern=".profile")
# print(bmi_files)

# Get list of INC PRS files
inc_files <- list.files(path="inc_prs_gctb/gctb/INC/INC_qc_n/score", pattern=".profile")
# print(inc_files)

##############

for(i in bmi_files){
	temp <- fread(paste0("inc_prs_gctb/gctb/BMI/BMI_qc_n/score/",i)) #%>% right_join(., sample.ids, by = c("FID" = "V1"))
	assign(gsub("\\..*","", i), temp)
}

bmi_full <- bmi_chr1

for(k in c(2:22)){
chri <- get(paste0("bmi_chr", k))
	names(chri)[6] <- paste0("SCORE_", k)
	bmi_full <- left_join(bmi_full, chri[, c(2,6)], by = "IID")
}
head(bmi_full)

##############

for(i in inc_files){
	temp <- fread(paste0("inc_prs_gctb/gctb/INC/INC_qc_n/score/",i)) #%>% right_join(., sample.ids, by = c("FID" = "V1"))
	assign(gsub("\\..*","", i), temp)
}

inc_full <- inc_chr1


# For subsequent chromosomes, add chromosomal SCORE to overall SCORE
for(k in c(2:22)){
	chri <- get(paste0("inc_chr", k))
n	ames(chri)[6] <- paste0("SCORE_", k)
	inc_full <- left_join(inc_full, chri[, c(2,6)], by = "IID")
}
head(inc_full)



build_betas <- data.frame(matrix(ncol=4, nrow=0))
names(build_betas) <- c("test_bmi_beta", "test_inc_beta", "build_bmi_beta", "build_inc_beta")

for(i in names(build)[4:97]){
	build_betas[i, "test_bmi_beta"] <- summary(lm(BMI ~ get(i), data=test ))$coef[2,1]
	build_betas[i, "test_inc_beta"] <- summary(lm(income ~ get(i), data=test ))$coef[2,1]
	build_betas[i, "build_bmi_beta"] <- summary(lm(BMI ~ get(i), data=build ))$coef[2,1]
	build_betas[i, "build_inc_beta"] <- summary(lm(income ~ get(i), data=build ))$coef[2,1]
}

build_se <- data.frame(matrix(ncol=4, nrow=0))
names(build_se) <- c("test_bmi_se", "test_inc_se", "build_bmi_se", "build_inc_se")

for(i in names(build)[4:97]){
	build_se[i, "test_bmi_se"] <- summary(lm(BMI ~ get(i), data=test ))$coef[2,2]
	build_se[i, "test_inc_se"] <- summary(lm(income ~ get(i), data=test ))$coef[2,2]
	build_se[i, "build_bmi_se"] <- summary(lm(BMI ~ get(i), data=build ))$coef[2,2]
	build_se[i, "build_inc_se"] <- summary(lm(income ~ get(i), data=build ))$coef[2,2]
}

##############



# Read complete UKBB Income dataset & subset height data
ukbb <- fread("../../ukbb_prs_complete_v2.txt")

all_snps <- fread("../../SNP_List/all_snps_v2.txt")
height <- subset(all_snps, Phenotype == "Height")
ukbb_ht <- select(ukbb, c(1, 17, 19,  which(names(ukbb) %in% height$SNP_1 | names(ukbb) %in% height$SNP_2)))
names(ukbb_ht)[4:3245] <- gsub("\\_.*", "", names(ukbb_ht)[4:3245])

# Read test set IDs
test_ids <- fread("inc_prs_gctb/gctb/test_subset_IDs.txt") 

# Create HT test subset
## use this when subsetting bootstrap samples within the loop
test <- ukbb_ht[ukbb_ht$ID %in% test_ids$V1, ]
# test <- ukbb[ukbb$ID %in% test_ids$V1, c(1, 18:19, 8549:8642)]
build <- ukbb_ht[!(ukbb_ht$ID %in% test_ids$V1), ]

# Get list of HEIGHT PRS files
ht_files <- list.files(path="inc_prs_gctb/gctb/HEIGHT/HEIGHT_qc_n/score/", pattern=".profile")
# print(ht_files)

# Get list of INC PRS files
inc_files <- list.files(path="inc_prs_gctb/gctb/INC/INC_qc_n/score", pattern=".profile")
# print(inc_files)

##############

for(i in ht_files){
	temp <- fread(paste0("inc_prs_gctb/gctb/HEIGHT/HEIGHT_qc_n/score/",i))
	assign(gsub("\\..*","", i), temp)
}

ht_full <- ht_chr1

### For subsequent chromosomes, add chromosomal SCORE to overall SCORE
for(k in c(2:22)){
	chri <- get(paste0("ht_chr", k))
	names(chri)[6] <- paste0("SCORE_", k)
	ht_full <- left_join(ht_full, chri[, c(2,6)], by = "IID")
}
head(ht_full)

##############

for(i in inc_files){
	temp <- fread(paste0("inc_prs_gctb/gctb/INC/INC_qc_n/score/",i))
	assign(gsub("\\..*","", i), temp)
}
	
inc_full <- inc_chr1

# For subsequent chromosomes, add chromosomal SCORE to overall SCORE
for(k in c(2:22)){
	chri <- get(paste0("inc_chr", k))
	names(chri)[6] <- paste0("SCORE_", k)
	inc_full <- left_join(inc_full, chri[, c(2,6)], by = "IID")
}
head(inc_full)

ht_snps_comp <- fread("../../SNP_List/Height_SNPs.txt") %>% 
select(., c("SNP", "CHR", "POS")) %>%
subset(., SNP %in% names(ukbb_ht)) %>%
mutate(., chr_cov = as.numeric(0), full_cov = as.numeric(0))

build_betas <- data.frame(matrix(ncol=4, nrow=0))
names(build_betas) <- c("test_ht_beta", "test_inc_beta", "build_ht_beta", "build_inc_beta")

for(i in names(build)[4:3245]){
	build_betas[i, "test_ht_beta"] <- summary(lm(height ~ get(i), data=test ))$coef[2,1]
	build_betas[i, "test_inc_beta"] <- summary(lm(income ~ get(i), data=test ))$coef[2,1]
	build_betas[i, "build_ht_beta"] <- summary(lm(height ~ get(i), data=build ))$coef[2,1]
	build_betas[i, "build_inc_beta"] <- summary(lm(income ~ get(i), data=build ))$coef[2,1]
}

build_se <- data.frame(matrix(ncol=4, nrow=0))
names(build_se) <- c("test_ht_se", "test_inc_se", "build_ht_se", "build_inc_se")

for(i in names(build)[4:3245]){
	build_se[i, "test_ht_se"] <- summary(lm(height ~ get(i), data=test ))$coef[2,2]
	build_se[i, "test_inc_se"] <- summary(lm(income ~ get(i), data=test ))$coef[2,2]
	build_se[i, "build_ht_se"] <- summary(lm(height ~ get(i), data=build ))$coef[2,2]
	build_se[i, "build_inc_se"] <- summary(lm(income ~ get(i), data=build ))$coef[2,2]
}


ht_snps <- ht_snps[, c("SNP", "CHR", "POS", "Tested_Allele", "Other_Allele", "test_ht_beta", "test_ht_se", "test_inc_beta", "test_inc_se", "build_ht_beta", "build_ht_se", "build_inc_beta", "build_inc_se")]
# 
# for(i in names(build)[4:3245]){
# build_betas[i, "test_ht_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(test[, ..i])), test[, "height"]))[2]
# build_betas[i, "test_inc_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(test[, ..i])), test[, "income"]))[2]
# build_betas[i, "build_ht_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(build[, ..i])), build[, "height"]))[2]
# build_betas[i, "build_inc_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(build[, ..i])), build[, "income"]))[2]
# rownames(build_betas)[i] <- names(build)[i]
# }