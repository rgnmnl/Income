##############################################################################
## Title: BMI PRS Bootstrap
## Version: 2
## Author: Regina Manansala
## Date Created: 06-November-2020
## Date Modified: 16-November-2020
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)

# Read complete UKBB Income dataset
ukbb <- fread("ukbb_prs_complete_v2.txt")

# Read test set IDs
test_ids <- fread("inc_prs_gctb/gctb/test_subset_IDs.txt") 

# Create BMI test subset
## use this when subsetting bootstrap samples within the loop
test <- ukbb[ukbb$ID %in% test_ids$V1, c(1, 18:19, 8549:8642)]

# Get list of BMI PRS files
bmi_files <- list.files(path="inc_prs_gctb/gctb/BMI/BMI_qc_n/score", pattern=".profile")
# print(bmi_files)

# Get list of INC PRS files
inc_files <- list.files(path="inc_prs_gctb/gctb/INC/INC_qc_n/score", pattern=".profile")
# print(inc_files)

#################################################################################################

## Initialize results matrix
bmi_betas_bs <- matrix(nrow=10000, ncol=94)
bmi_betas_adj <- matrix(nrow=10000, ncol=94)

#################################################################################################

## Initialize results matrix
inc_betas_bs <- matrix(nrow=10000, ncol=94)
inc_betas_adj <- matrix(nrow=10000, ncol=94)

################# 

start <- Sys.time()
print(paste0("Start: ", start))
for(b in 1:100){

	## Sample
	sample.ids <- sample(test_ids$V1, size=length(test_ids$V1), replace=TRUE) %>% as.data.frame() %>% rename(., "V1" = `.`)

	################# BMI - Calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score

	## Get BMI SNPs from sBayesR PRS results
	bmi_snps_comp <- fread("SNP_List/BMI_SNPs.txt") %>% 
		select(., -c("Allele1", "Allele2")) %>%
		mutate(., chr_cov = as.numeric(0), full_cov = as.numeric(0))
	# head(bmi_snps_comp)

	## Read BMI PRS files and subset to only include test samples
	for(i in bmi_files){
		temp <- fread(paste0("inc_prs_gctb/gctb/BMI/BMI_qc_n/score/",i)) %>% right_join(., sample.ids, by = c("FID" = "V1"))
		assign(gsub("\\..*","", i), temp)
	}

	## Calculate Overall PRS score
	bmi_full <- bmi_chr1
	
	### For subsequent chromosomes, add chromosomal SCORE to overall SCORE
	for(k in c(2:22)){
		chri <- get(paste0("bmi_chr", k))
		bmi_full$SCORE <- bmi_full$SCORE + chri$SCORE
	}
	head(bmi_full)

	## calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score
	for(i in 1:nrow(bmi_snps_comp)){
		snp <- test[[bmi_snps_comp$SNP[i]]]
		chr <- bmi_snps_comp$CHR[i]
	
		chr_prs <- get(paste0("bmi_chr", chr))$SCORE
	
		bmi_snps_comp[i, 5] <- cov(snp, chr_prs)
		bmi_snps_comp[i, 6] <- cov(snp, bmi_full$SCORE)
	}
	# head(bmi_snps_comp)

	# Calculate BMI adjustments
	bmi_adjustments <- bmi_snps_comp$chr_cov/bmi_snps_comp$full_cov

	bmi <- left_join(sample.ids, test, by = c("V1" = "ID"))

	build_betas <- 	data.frame(matrix(ncol=1, nrow=0))
	names(build_betas) <- c("bmi_beta")

	for(i in names(bmi)[4:97]){
		build_betas[i, "bmi_beta"] <- summary(lm(BMI ~ get(i), data=bmi ))$coef[2,1]
	}

	build_betas$new_bmi_beta <- ifelse(build_betas$bmi_beta < 0, (-1)*build_betas$bmi_beta, build_betas$bmi_beta)
	bmi_betas_bs[b, ] <- t(build_betas$bmi_beta)

	build_betas$adjusted.bmi.beta <- build_betas$new_bmi_beta*bmi_adjustments
	bmi_betas_adj[b, ] <- t(build_betas$adjusted.bmi.beta)

	################# INCOME - Calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score

	## Get INC SNPs from sBayesR PRS results
	inc_snps_comp <- fread("SNP_List/BMI_SNPs.txt") %>% 
		select(., -c("Allele1", "Allele2")) %>%
		mutate(., chr_cov = as.numeric(0), full_cov = as.numeric(0))
	# head(inc_snps_comp)

	## Read INC PRS files and subset to only include test samples
	for(i in inc_files){
		temp <- fread(paste0("inc_prs_gctb/gctb/INC/INC_qc_n/score/",i)) %>% right_join(., sample.ids, by = c("FID" = "V1"))
		assign(gsub("\\..*","", i), temp)
	}

	## Calculate Overall PRS score
	inc_full <- inc_chr1
	
	# For subsequent chromosomes, add chromosomal SCORE to overall SCORE
	for(k in c(2:22)){
		chri <- get(paste0("inc_chr", k))
		inc_full$SCORE <- inc_full$SCORE + chri$SCORE
	}

	## calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score
	for(i in 1:nrow(inc_snps_comp)){
		snp <- test[[inc_snps_comp$SNP[i]]]
		chr <- inc_snps_comp$CHR[i]
	
		chr_prs <- get(paste0("inc_chr", chr))$SCORE
	
		inc_snps_comp[i, 5] <- cov(snp, chr_prs)
		inc_snps_comp[i, 6] <- cov(snp, inc_full$SCORE)
	}
	# head(inc_snps_comp)

	inc_adjustments <- inc_snps_comp$chr_cov/inc_snps_comp$full_cov

	################################## 
	
	inc <- left_join(sample.ids, test, by = c("V1" = "ID"))
	
	inc_build_betas <- 	data.frame(matrix(ncol=1, nrow=0))
	names(inc_build_betas) <- c("inc_beta")

	for(i in names(inc)[4:97]){
		inc_build_betas[i, "inc_beta"] <- summary(lm(income ~ get(i), data=inc ))$coef[2,1]
	}
	
	inc_build_betas$new_income_beta <- ifelse(build_betas$bmi_beta < 0, (-1)*inc_build_betas$inc_beta, inc_build_betas$inc_beta)
	inc_betas_bs[b, ] <- t(inc_build_betas$inc_beta)
	
	inc_build_betas$adjusted.inc.beta <- inc_build_betas$new_income_beta*inc_adjustments
	inc_betas_adj[b, ] <- t(inc_build_betas$adjusted.inc.beta)
}
end <- Sys.time()
print(paste0("End: ", end))

print(end-start)

colnames(bmi_betas_bs) <- bmi_snps_comp$SNP
colnames(bmi_betas_adj) <- bmi_snps_comp$SNP
colnames(inc_betas_bs) <- bmi_snps_comp$SNP
colnames(inc_betas_adj) <- bmi_snps_comp$SNP

write.table(bmi_betas_bs, "BOOTSTRAP/BMI_betas_bs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(bmi_betas_adj, "BOOTSTRAP/BMI_betas_adj.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

write.table(inc_betas_bs, "BOOTSTRAP/INC_betas_bs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(inc_betas_adj, "BOOTSTRAP/INC_betas_adj.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)




























