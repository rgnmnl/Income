##############################################################################
## Title: sBayesR Variance, Covariance, Correlation function
## Author: Regina Manansala
## Date Created: 14-October-2020
## Date Modified: 04-November-2020
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

sbr_var <- function(phenotype, ids = FALSE){

	# Specify the phenotype
	i <- phenotype
	
	# Specify which set of sample IDs to use (pre-determined test group vs. 157k sampled IDs from test group)
	if(ids == FALSE){
		IDs <- fread("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/test_subset_IDs.txt") #%>% `[[`("V1") %>% as.list() %>% unlist()
	} else {
		IDs <- sample.ids
	}
	
	# Get prefix of score file for specified phenotype
	file <- gsub("\\_.*","",list.files(path=paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/",i,"/", i, "_qc_n/score"), pattern= "chr1.profile"))

	# Initialize data frame to store individual chromosome variances
	variances <- data.frame(matrix(ncol=2, nrow=0))
	names(variances) <- c("type", "var")

	# For each chromosome, read PLINK score file and subset based on specified sample IDs
	# Calculate variance for each chromosome and store in variances data frame
	# Output subsetted PLINK chromosome score files
	for(j in 1:22){
# 	for(j in c(1:5,7:22)){
		temp <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/",i,"/", i, "_qc_n/score/", file, "_chr", j, ".profile")) %>% right_join(., IDs, by = c("FID" = "V1")) #%>% subset(., FID %in% IDs)
		variances[j, "type"] <- paste0("chr", j)
		variances[j, "var"] <- var(temp$SCORE)
		assign(paste0("chr", j), temp)
	}

	# Initialize overall PRS data frame with chr 1
	full <- chr1
	
	# For subsequent chromosomes, add chromosomal SCORE to overall SCORE
	for(k in seq(2, 22)){
# 	for(k in c(2:5,7:22)){
		chri <- get(paste0("chr", k))
		full$SCORE <- full$SCORE + chri$SCORE
	}

	# Calculate variance of PRS:
	# (Variance of sum of all scores - sum of chromosomal variances) / Variance of sum of all scores
	(var(full$SCORE) - sum(variances$var, na.rm=TRUE)) / var(full$SCORE)
}

# sbr_var("BMI")
# sbr_var("INC")
# sbr_var("HEIGHT")
# sbr_var("WBC")
# sbr_var("RBC")
# sbr_var("EDU_num")

sbr_cov <- function(pheno1, pheno2, ids = FALSE){
	
	# Specify the phenotypes
	a <- pheno1
	b <- pheno2
	
	# Specify which set of sample IDs to use (pre-determined test group vs. 157k sampled IDs from test group)
	if(ids == FALSE){
		IDs <- fread("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/test_subset_IDs.txt") #%>% `[[`("V1") %>% as.list() %>% unlist()
	} else {
		IDs <- sample.ids
	}
	
	# Get prefixes of score files for specified phenotypes
	file_a <- gsub("\\_.*","",list.files(path=paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/",a,"/", a, "_qc_n/score"), pattern= "chr1.profile"))
	file_b <- gsub("\\_.*","",list.files(path=paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/",b,"/", b, "_qc_n/score"), pattern= "chr1.profile"))

	# Initialize data frame to store chromosomal covariances
	covar <- data.frame(matrix(ncol=2, nrow=0))
	names(covar) <- c("type", "covar")

	# For each chromosome, read PLINK score files and subset based on specified sample IDs
	# Calculate covariance between the 2 phenotypes for each chromosome and store in covariances data frame
	# Output subsetted PLINK chromosome score files
	for(j in 1:22){
# 	for(j in c(1:5,7:22)){
		temp_a <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/", a,"/", a, "_qc_n/score/", file_a, "_chr", j, ".profile")) %>% right_join(., IDs, by = c("FID" = "V1")) #%>% subset(., FID %in% IDs)
		temp_b <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/", b,"/", b, "_qc_n/score/", file_b, "_chr", j, ".profile")) %>% right_join(., IDs, by = c("FID" = "V1")) #%>% subset(., FID %in% IDs)
		
		covar[j, "type"] <- paste0("chr", j)
		covar[j, "covar"] <- cov(temp_a$SCORE, temp_b$SCORE)
		
		assign(paste0(file_a, "_chr", j), temp_a)
		assign(paste0(file_b, "_chr", j), temp_b)
	}
	
	# Initialize overall PRS data frame with chr 1 for each phenotype
	full_a <- get(paste0(file_a, "_chr1"))
	full_b <- get(paste0(file_b, "_chr1"))
	
	# For subsequent chromosomes, add chromosomal SCORE to overall SCORE
	for(k in seq(2, 22)){
# 	for(k in c(2:5,7:22)){
		chra <- get(paste0(file_a, "_chr", k))
		full_a$SCORE <- full_a$SCORE + chra$SCORE

		chrb <- get(paste0(file_b, "_chr", k))
		full_b$SCORE <- full_b$SCORE + chrb$SCORE
	}
	# Calculate covariance of PRS:
	# (Covariance of sum of all scores - sum of chromosomal covariances) / Covariance of sum of all scores
	(cov(full_a$SCORE, full_b$SCORE) - sum(covar$covar, na.rm=TRUE)) / cov(full_a$SCORE, full_b$SCORE)
}

# sbr_cov("EDU", "INC")
# sbr_cov("HEIGHT", "INC")
# sbr_cov("BMI", "INC")
# sbr_cov("RBC", "WBC")
# sbr_cov("RBC", "INC")
# sbr_cov("WBC", "INC")

sbr_corr <- function(pheno1, pheno2, ids = FALSE){

	# Specify the phenotypes
	a <- pheno1
	b <- pheno2
	
	# Specify which set of sample IDs to use (pre-determined test group vs. 157k sampled IDs from test group)
	if(ids == FALSE){
		IDs <- fread("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/test_subset_IDs.txt") #%>% `[[`("V1") %>% as.list() %>% unlist()
	} else {
		IDs <- sample.ids
	}
	
	# Get prefixes of score files for specified phenotypes
	file_a <- gsub("\\_.*","",list.files(path=paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/",a,"/", a, "_qc_n/score"), pattern= "chr1.profile"))
	file_b <- gsub("\\_.*","",list.files(path=paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/",b,"/", b, "_qc_n/score"), pattern= "chr1.profile"))

	# For each chromosome, read PLINK score files and subset based on specified sample IDs
	# Output subsetted PLINK chromosome score files
	for(j in 1:22){
		temp_a <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/", a,"/", a, "_qc_n/score/", file_a, "_chr", j, ".profile")) %>% right_join(., IDs, by = c("FID" = "V1"))
		temp_b <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs_gctb/gctb/", b,"/", b, "_qc_n/score/", file_b, "_chr", j, ".profile")) %>% right_join(., IDs, by = c("FID" = "V1"))

		assign(paste0(file_a, "_chr", j), temp_a)
		assign(paste0(file_b, "_chr", j), temp_b)
	}
	
	# Initialize overall PRS data frame with chr 1 for each phenotype
	full_a <- get(paste0(file_a, "_chr1"))
	full_b <- get(paste0(file_b, "_chr1"))
	
	# For subsequent chromosomes, add chromosomal SCORE to overall SCORE
	for(k in seq(2, 22)){
		chra <- get(paste0(file_a, "_chr", k))
		full_a$SCORE <- full_a$SCORE + chra$SCORE

		chrb <- get(paste0(file_b, "_chr", k))
		full_b$SCORE <- full_b$SCORE + chrb$SCORE
	}
	
	# Calculate correlation of PRS:
	cor(full_a$SCORE, full_b$SCORE)

}

# sbr_corr("EDU", "INC")
# sbr_corr("HEIGHT", "INC")
# sbr_corr("BMI", "INC")
# sbr_corr("RBC", "WBC")
# sbr_corr("RBC", "INC")
# sbr_corr("WBC", "INC")







