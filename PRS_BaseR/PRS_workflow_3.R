##############################################################################
## Title: PRS Workflow Step 3
## Purpose: PLINK score analysis
## Author: Regina Manansala
## Date Created: 7-August-2019
## Date Modified: 22-June-2020
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

var <- "dep"

build <- fread(paste0("../PLINK_DEP/", var, "_pheno_build.txt"))
test <- fread(paste0("../PLINK_DEP/", var, "_pheno_test.txt"))


build <- fread(paste0(var, "_pheno_build.txt"))
test <- fread(paste0(var, "_pheno_test.txt"))

score_files <- list.files(pattern="chr([0-9]{1,2})_test.txt")
# score_files <- gsub("\\..*","", score_files)

for(i in score_files){
	assign(gsub("\\..*","", i), fread(i))
}