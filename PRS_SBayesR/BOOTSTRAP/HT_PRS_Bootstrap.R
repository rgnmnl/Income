##############################################################################
## Title: Height PRS Bootstrap
## Author: Regina Manansala
## Date Created: 04-December-2020
## Date Modified: 05-February-2021
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)

# Read complete UKBB Income dataset & subset height data
ukbb <- fread("ukbb_prs_complete_v2.txt")

all_snps <- fread("SNP_List/all_snps_v2.txt")
height <- subset(all_snps, Phenotype == "Height")
ukbb_ht <- select(ukbb, c(1, 17, 19,  which(names(ukbb) %in% height$SNP_1 | names(ukbb) %in% height$SNP_2)))
names(ukbb_ht)[4:3245] <- gsub("\\_.*", "", names(ukbb_ht)[4:3245])

# Read test set IDs
test_ids <- fread("inc_prs_gctb/gctb/test_subset_IDs.txt") 

# Create HT test subset
## use this when subsetting bootstrap samples within the loop
test <- ukbb_ht[ukbb_ht$ID %in% test_ids$V1, ]
# test <- ukbb[ukbb$ID %in% test_ids$V1, c(1, 18:19, 8549:8642)]

# Get list of HEIGHT PRS files
ht_files <- list.files(path="inc_prs_gctb/gctb/HEIGHT/HEIGHT_qc_n/score/", pattern=".profile")
# print(ht_files)

# Get list of INC PRS files
inc_files <- list.files(path="inc_prs_gctb/gctb/INC/INC_qc_n/score", pattern=".profile")
# print(inc_files)

#################################################################################################

ht_betas_bs <- matrix(nrow=10000, ncol=3242)
ht_betas_adj <- matrix(nrow=10000, ncol=3242)

#################################################################################################

inc_betas_bs <- matrix(nrow=10000, ncol=3242)
inc_betas_adj <- matrix(nrow=10000, ncol=3242)

#################################################################################################

# Betas without sampling

## Get HEIGHT SNPs from sBayesR PRS results
ht_snps_comp <- fread("SNP_List/Height_SNPs.txt") %>% 
	select(., c("SNP", "CHR", "POS")) %>%
	subset(., SNP %in% names(ukbb_ht)) %>%
	mutate(., chr_cov = as.numeric(0), full_cov = as.numeric(0))

## Read HEIGHT PRS files and subset to only include test samples
for(i in ht_files){
	temp <- fread(paste0("inc_prs_gctb/gctb/HEIGHT/HEIGHT_qc_n/score/",i)) %>% subset(., FID %in% test_ids$V1)
	assign(gsub("\\..*","", i), temp)
}

## Calculate Overall PRS score
ht_full <- ht_chr1
	
### For subsequent chromosomes, add chromosomal SCORE to overall SCORE
for(k in c(2:22)){
	chri <- get(paste0("ht_chr", k))
	ht_full$SCORE <- ht_full$SCORE + chri$SCORE
}
head(ht_full)

## calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score
for(i in 1:nrow(ht_snps_comp)){
	snp <- test[[ht_snps_comp$SNP[i]]]
	chr <- ht_snps_comp$CHR[i]
	
	chr_prs <- get(paste0("ht_chr", chr))$SCORE
	
	ht_snps_comp[i, 4] <- cov(snp, chr_prs)
	ht_snps_comp[i, 5] <- cov(snp, ht_full$SCORE)
} ###### *********** ######

ht_adjustments <- ht_snps_comp$chr_cov/ht_snps_comp$full_cov

test <- as.data.frame(test)

og_build_betas <- 	data.frame(matrix(ncol=1, nrow=0))
names(og_build_betas) <- c("ht_beta")

for(i in 1:3242){
	og_build_betas[i, "ht_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(test[, i + 3])), test[, "height"]))[2]	
	rownames(og_build_betas)[i] <- names(test)[i + 3]
}

# for(i in names(test)[4:3245]){
# 	og_build_betas[i, "ht_beta_reg_lm"] <- summary(lm(height ~ get(i), data=test))$coef[2,1]
# }

og_build_betas$new_ht_beta <- ifelse(og_build_betas$ht_beta < 0, (-1)*og_build_betas$ht_beta, og_build_betas$ht_beta)
og_build_betas$adjusted.ht.beta <- og_build_betas$new_ht_beta*ht_adjustments

write.table(og_build_betas, "BOOTSTRAP/HEIGHT/HT_betas_og.txt", sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#################################################################################################


## Get INC SNPs from sBayesR PRS results
inc_snps_comp <- fread("SNP_List/Height_SNPs.txt") %>% 
	select(., c("SNP", "CHR", "POS")) %>%
	subset(., SNP %in% names(ukbb_ht)) %>%
	mutate(., chr_cov = as.numeric(0), full_cov = as.numeric(0))

## Read INC PRS files and subset to only include test samples
for(i in inc_files){
	temp <- fread(paste0("inc_prs_gctb/gctb/INC/INC_qc_n/score/",i)) %>% subset(., FID %in% test_ids$V1)
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
	
	inc_snps_comp[i, 4] <- cov(snp, chr_prs)
	inc_snps_comp[i, 5] <- cov(snp, inc_full$SCORE)
}
	# head(inc_snps_comp)

inc_adjustments <- inc_snps_comp$chr_cov/inc_snps_comp$full_cov

og_inc_build_betas <- 	data.frame(matrix(ncol=1, nrow=0))
names(og_inc_build_betas) <- c("inc_beta")

for(i in 1:3242){
	og_inc_build_betas[i, "inc_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(test[, i + 3])), test[, "income"]))[2]	
	rownames(og_inc_build_betas)[i] <- names(test)[i + 3]
}

# for(i in names(test)[4:3245]){
# 	og_inc_build_betas[i, "inc_beta_reg_lm"] <- summary(lm(income ~ get(i), data=test))$coef[2,1]
# }

og_inc_build_betas$new_inc_beta <- ifelse(og_inc_build_betas$inc_beta < 0, (-1)*og_inc_build_betas$inc_beta, og_inc_build_betas$inc_beta)
og_inc_build_betas$adjusted.inc.beta <- og_inc_build_betas$new_inc_beta*inc_adjustments

write.table(og_inc_build_betas, "BOOTSTRAP/HEIGHT/INC_HT_betas_og.txt", sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)


#################

start <- Sys.time()
print(paste0("Start: ", start))
for(b in 1:10){
	sample.ids <- sample(test_ids$V1, size=length(test_ids$V1), replace=TRUE) %>% as.data.frame() %>% rename(., "V1" = `.`)

	################# HEIGHT - Calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score

	## Get HEIGHT SNPs from sBayesR PRS results
	ht_snps_comp <- fread("SNP_List/Height_SNPs.txt") %>% 
		select(., c("SNP", "CHR", "POS")) %>%
		subset(., SNP %in% names(ukbb_ht)) %>%
		mutate(., chr_cov = as.numeric(0), full_cov = as.numeric(0))
	# head(bmi_snps_comp)
	
	## Read HEIGHT PRS files and subset to only include test samples
	for(i in ht_files){
		temp <- fread(paste0("inc_prs_gctb/gctb/HEIGHT/HEIGHT_qc_n/score/",i)) %>% right_join(., sample.ids, by = c("FID" = "V1"))
		assign(gsub("\\..*","", i), temp)
	}

	## Calculate Overall PRS score
	ht_full <- ht_chr1
	
	### For subsequent chromosomes, add chromosomal SCORE to overall SCORE
	for(k in c(2:22)){
		chri <- get(paste0("ht_chr", k))
		ht_full$SCORE <- ht_full$SCORE + chri$SCORE
	}
	head(ht_full)

	## calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score
# 	start <- Sys.time()
	for(i in 1:nrow(ht_snps_comp)){
		snp <- test[[ht_snps_comp$SNP[i]]]
		chr <- ht_snps_comp$CHR[i]
	
		chr_prs <- get(paste0("ht_chr", chr))$SCORE
	
		ht_snps_comp[i, 4] <- cov(snp, chr_prs)
		ht_snps_comp[i, 5] <- cov(snp, ht_full$SCORE)
	} ###### *********** ######

	# Calculate HT adjustments
	ht_adjustments <- ht_snps_comp$chr_cov/ht_snps_comp$full_cov

	ht <- left_join(sample.ids, test, by = c("V1" = "ID"))

	build_betas <- 	data.frame(matrix(ncol=1, nrow=0))
	names(build_betas) <- c("ht_beta")

	###### *********** ######
# 	for(i in names(ht)[4:3245]){
# 		build_betas[i, "ht_beta"] <- summary(lm(height ~ get(i), data=ht))$coef[2,1]
# 	}

	for(i in 1:3242){
		build_betas[i, "ht_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(ht[, i + 3])), ht[, "height"]))[2]	
		rownames(build_betas)[i] <- names(ht)[i + 3]
	}

	###### *********** ######

	build_betas$new_ht_beta <- ifelse(build_betas$ht_beta < 0, (-1)*build_betas$ht_beta, build_betas$ht_beta)
	ht_betas_bs[b, ] <- t(build_betas$ht_beta)

	build_betas$adjusted.ht.beta <- build_betas$new_ht_beta*ht_adjustments
	ht_betas_adj[b, ] <- t(build_betas$adjusted.ht.beta)

	################# INCOME - Calculate covariance between SNP genotype and chromosomal score && SNP genotype and overall score

	## Get INC SNPs from sBayesR PRS results
	inc_snps_comp <- fread("SNP_List/Height_SNPs.txt") %>% 
		select(., c("SNP", "CHR", "POS")) %>%
		subset(., SNP %in% names(ukbb_ht)) %>%
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
	
		inc_snps_comp[i, 4] <- cov(snp, chr_prs)
		inc_snps_comp[i, 5] <- cov(snp, inc_full$SCORE)
	}
	# head(inc_snps_comp)

	inc_adjustments <- inc_snps_comp$chr_cov/inc_snps_comp$full_cov

	################################## 
	
	inc <- left_join(sample.ids, test, by = c("V1" = "ID"))
	
	inc_build_betas <- 	data.frame(matrix(ncol=1, nrow=0))
	names(inc_build_betas) <- c("inc_beta")

# 	for(i in names(inc)[4:3245]){
# 		inc_build_betas[i, "inc_beta"] <- summary(lm(income ~ get(i), data=inc))$coef[2,1]
# 	}
	
	for(i in 1:3242){
		inc_build_betas[i, "inc_beta"] <- coef(lm.fit(x = cbind(1,as.matrix(inc[, i + 3])), inc[, "income"]))[2]	
		rownames(inc_build_betas)[i] <- names(inc)[i + 3]
	}
	
	inc_build_betas$new_income_beta <- ifelse(build_betas$ht_beta < 0, (-1)*inc_build_betas$inc_beta, inc_build_betas$inc_beta)
	inc_betas_bs[b, ] <- t(inc_build_betas$inc_beta)
	
	inc_build_betas$adjusted.inc.beta <- inc_build_betas$new_income_beta*inc_adjustments
	inc_betas_adj[b, ] <- t(inc_build_betas$adjusted.inc.beta)
}
end <- Sys.time()
print(paste0("End: ", end))

print(end-start)

colnames(ht_betas_bs) <- ht_snps_comp$SNP
colnames(ht_betas_adj) <- ht_snps_comp$SNP
colnames(inc_betas_bs) <- ht_snps_comp$SNP
colnames(inc_betas_adj) <- ht_snps_comp$SNP

write.table(ht_betas_bs, "BOOTSTRAP/HT_betas_bs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(ht_betas_adj, "BOOTSTRAP/HT_betas_adj.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

write.table(inc_betas_bs, "BOOTSTRAP/INC_HT_betas_bs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(inc_betas_adj, "BOOTSTRAP/INC_HT_betas_adj.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)













