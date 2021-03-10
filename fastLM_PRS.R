##############################################################################
## Title: Polygenic Risk Score Calculation using fastLM
## Author: Regina Manansala
## Date Created: 01-April-2020
## Date Modified: 04-May-2020
##############################################################################

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)
library(tibble)
library(stringr)

## Import BED files (LD pruned and MAF filtered in plink)
chr1 <- BEDMatrix("/mafsub_chr1.bed")
chr2 <- BEDMatrix("/mafsub_chr2.bed")
chr3 <- BEDMatrix("/mafsub_chr3.bed")
chr4 <- BEDMatrix("/mafsub_chr4.bed")
chr5 <- BEDMatrix("/mafsub_chr5.bed")
chr6 <- BEDMatrix("/mafsub_chr6.bed")
chr7 <- BEDMatrix("/mafsub_chr7.bed")
chr8 <- BEDMatrix("/mafsub_chr8.bed")
chr9 <- BEDMatrix("/mafsub_chr9.bed")
chr10 <- BEDMatrix("/mafsub_chr10.bed")
chr11 <- BEDMatrix("/mafsub_chr11.bed")
chr12 <- BEDMatrix("/mafsub_chr12.bed")
chr13 <- BEDMatrix("/mafsub_chr13.bed")
chr14 <- BEDMatrix("/mafsub_chr14.bed")
chr15 <- BEDMatrix("/mafsub_chr15.bed")
chr16 <- BEDMatrix("/mafsub_chr16.bed")
chr17 <- BEDMatrix("/mafsub_chr17.bed")
chr18 <- BEDMatrix("/mafsub_chr18.bed")
chr19 <- BEDMatrix("/mafsub_chr19.bed")
chr20 <- BEDMatrix("/mafsub_chr20.bed")
chr21 <- BEDMatrix("/mafsub_chr21.bed")

## Combine bed files into matrix
wg <- ColumnLinkedMatrix(chr1, chr2, chr3, chr4, chr5, chr6, chr7,
						chr8, chr9, chr10, chr11, chr12, chr13, chr14,
						chr15, chr16, chr17, chr18, chr19, chr20, chr21)

## Impute missing genotype date with mean		
# imp <- wg %>% as.matrix()
# imp[is.na(wg)] <- mean(wg, na.rm=TRUE)

## Import Phenotype file
# dat <- fread("ld_ukbb_prs_complete.txt", header=TRUE)
# pheno <- dat[, 1:21]
# rownames(pheno) <- paste(pheno$ID, pheno$ID, sep = "_")

pheno <- fread("build_pheno.txt")
rownames(pheno) <- paste(pheno$IID, pheno$IID, sep = "_")

## Clear BED files
# rm(dat)
rm(list=ls(pattern="chr"))
gc()

## Function for running GWAS
## Adjust for Age and PCs
## 	1. Regress ht vs Age + PCs

## Regression using BGData package
# wg_sub <- wg[, 1:100000]
# DATA <- BGData(geno = wg_sub, pheno = pheno)
# GWAS <- GWAS(formula = height ~ 1 + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = DATA)

# DATA <- BGData(geno = wg, pheno = pheno)
# GWAS <- GWAS(formula = height ~ 1, data = DATA)

## Regression using base R
mod <- lm(height ~ Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno)

## 	2. Get residuals

y <- resid(mod)

## 	3. These residuals will be the y-variable that's input to the fast_lm function
## 	4. X is a matrix with 2 columns. First column is all 1s (intercept). 2nd column is genotypes.

## Initialize fastLM function
fast_lm <- function(y,X){
    n<-nrow(X);
    A <- t(X)%*%X;
    Ai <- solve(A);
    mat1 <- Ai%*%t(X)
    beta <- mat1%*%y
    y.hat <- X%*%beta;
    resid <- y-y.hat;
    covar <- Ai*as.numeric(t(resid)%*%resid)/(n-ncol(X));
    res_list <- vector('list',6);
    names(res_list)<-c('mat1','resid','se','y.hat','z', 'beta');
    res_list$mat1<-mat1;
    res_list$resid<-resid;
    res_list$se<-diag(Ai);
    res_list$y.hat <- y.hat
    res_list$z <- beta/sqrt(diag(covar))
    res_list$beta <- beta
    res_list$p <- pnorm(abs(res_list$z), lower.tail=FALSE)
    return(res_list);
}

## 
start <- proc.time()
for(i in 1:10000){
	wg_sub <- wg[,i] %>% as.matrix() %>% as.data.frame() %>% 
		rownames_to_column("ID") %>%
		mutate(., V1 = as.numeric(V1)) %>%
	# 	mutate(., ID = gsub("_.*", "", ID), V1 = as.numeric(V1)) %>%
		add_column(., Intercept = 1, .before = "V1") %>%
		column_to_rownames("ID")
	wg_sub$V1 <- ifelse(is.na(wg_sub$V1), mean(wg_sub$V1, na.rm=TRUE), wg_sub$V1)
	colnames(wg_sub)[2] <- colnames(wg)[i]
	
	X <- as.matrix(wg_sub[rownames(wg_sub) %in% names(y),])
	
	foo <- fast_lm(y = y, X= X)
	assign(paste0("out_", i), foo)
	print(i)
}
end <- proc.time()

save(list = ls(pattern = "out"), file = "Height_fastlm_out.RData")