##############################################################################
## Title: UKBB MI GCTA
## Note: Run MI gwas using GCTA
## Author: Regina Manansala
## Date Created: 14-February-2021
## Date Modified: 18-February-2021
##############################################################################

###############################################
######## Run MI sBr using LDL SNP info ########
###############################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

## Read MI SAIGE results
files <- list.files(path="gctb/MI/SAIGE/", pattern=".SAIGE.bgen.txt")

full <- data.table()
for(i in files){
	temp <- fread(paste0("gctb/MI/SAIGE/", i))
	full <- rbind(full, temp)
}

## Read LDL GCTA filtered summary stats
qc <- c("", "_pval_0.1", "_pval_0.2", "_pval_0.3", "_pval_0.4", "_pval_0.5", "_pval_0.6", "_pval_0.7", "_pval_0.8")

for(i in qc){
	temp <- fread(paste0("gctb/LDL/LDL_qc_n", i, "/ldl_qc_n", i, ".ma")) %>% 
		left_join(., full[, c("rsid", "Allele1", "Allele2", "BETA", "SE", "p.value")], by = c("SNP"="rsid", "A1" = "Allele1", "A2"= "Allele2")) %>%
		select(., c("SNP", "A1", "A2", "AF1", "BETA.y", "SE.y", "p.value", "N")) %>%
		rename(., BETA = BETA.y, SE = SE.y, P = p.value)
	write.table(temp, paste0("mi_ldl_qc_n", i, ".txt"), sep = "\t", row.names = FALSE, col.names=TRUE, quote = FALSE)
# 	assign(paste0("mi_ldl_qc_n", i), temp)
}

## Read MI GCTA filtered summary stats
qc <- c("", "_pval_0.1", "_pval_0.2", "_pval_0.3", "_pval_0.4", "_pval_0.5", "_pval_0.6", "_pval_0.7", "_pval_0.8")

for(i in qc){
	temp <- fread(paste0("gctb/MI/MI_qc_n", i, "/mi_qc_n", i, ".ma")) %>%
		rename(., AF1 = MAF)
	write.table(temp[, c("SNP", "A1", "A2", "AF1", "BETA", "SE", "P", "N")], paste0("gctb/MI/MI_qc_n", i, "/mi_qc_n", i, ".ma"), sep = "\t", row.names = FALSE, col.names=TRUE, quote = FALSE)
# 	assign(paste0("mi_ldl_qc_n", i), temp)
}

# qc <- c("", "_pval_0.1", "_pval_0.2", "_pval_0.3", "_pval_0.4", "_pval_0.5", "_pval_0.6", "_pval_0.7", "_pval_0.8")
# 
# for(i in qc){
# 	temp <- fread(paste0("gctb/MI/MI_qc_n", i, "/mi_qc_n", i, ".ma"))
# 	temp <- temp[,c("SNP", "A1", "A2", "AF1", "BETA", "SE", "P", "N")]
# 	assign(paste0("mi_qc_n", i), temp)
# # 	write.table(temp, paste0("../MI_qc_n", i, "/mi_qc_n", i, ".ma"), col.names = T, row.names = F, sep = " ", quote = F)
# }

###############################################################
######### Step 3 - Make params file for plink scoring #########
###############################################################

library(data.table)
library(dplyr)

## Read LDL/MI sBayesR results
x <- read.table("TEST_qc_n/mi_ldl_qc_n.snpRes", header = T)
x.2 <- x[,c("Name", "A1", "A1Effect")]
write.table(x.2, "TEST_qc_n/mi_ldl_qc_n.params", row.names = F, sep = " ", quote = F)

for(i in 1:8){
	x <- read.table(paste0("TEST_qc_n_pval_0.", i, "/mi_ldl_qc_n_pval_0.", i, ".snpRes"), header = T)
	x.2 <- x[,c("Name", "A1", "A1Effect")]
	write.table(x.2, paste0("TEST_qc_n_pval_0.", i, "/mi_ldl_qc_n_pval_0.", i, ".params"), row.names = F, sep = " ", quote = F)
}

###############################################
######### Step 5 - Calculate Pred R^2 #########
###############################################

library(data.table)
library(dplyr)

ukb <- read.table("gctb/MI/SAIGE/ukbb_prs_pheno_mi_test.txt", header = TRUE)

resLocke <- data.frame(matrix(0, nrow = 9, ncol = 2))
qc <- c("", "_pval_0.1", "_pval_0.2", "_pval_0.3", "_pval_0.4", "_pval_0.5", "_pval_0.6", "_pval_0.7", "_pval_0.8")

for(k in 1:length(qc)){
	lky.chr1 <- read.table(paste0("TEST_qc_n", qc[k], "/score/mi_ldl_chr1.profile"), header = T)
	lky.all  <- lky.chr1
	for (i in seq(2, 22)){
# 	for(i in c(2:5, 7:22)){
		print(paste("Doing", i))
		lky.chri <- read.table(paste0("TEST_qc_n", qc[k], "/score/mi_ldl_chr", i, ".profile"), header = T)
		lky.all$SCORE <- lky.all$SCORE + lky.chri$SCORE
  	}
  	lky.ukb <- inner_join(ukb, lky.all, by = c("ID" = "FID"))
  	resLocke[k, 1] <- paste0("TEST_qc_n", qc[k])
	resLocke[k, 2] <- summary(lm(lky.ukb$mi_flag ~ lky.ukb$SCORE))$r.squared
	print(resLocke)
}
colnames(resLocke) <- c("QC", "Pred_R2")

write.table(resLocke, "mi_ldl_test_results.txt", col.names = T, row.names = F, quote = F, sep = " ")


#####################################################################
############## Get plink scores from CAD Summary Stats ##############
#####################################################################

##################################################
############## dbSNP Data Wrangling ##############
##################################################

awk '{gsub(/NC_000001.10/, "1"); gsub(/NC_000002.11/, "2"); gsub(/NC_000003.11/, "3"); gsub(/NC_000004.11/, "4"); gsub(/NC_000005.9/, "5"); gsub(/NC_000006.11/, "6"); gsub(/NC_000007.13/, "7"); gsub(/NC_000008.10/, "8"); gsub(/NC_000009.11/, "9"); gsub(/NC_000010.10/, "10"); gsub(/NC_000011.9/, "11"); gsub(/NC_000012.11/, "12"); gsub(/NC_000013.10/, "13"); gsub(/NC_000014.8/, "14"); gsub(/NC_000015.9/, "15"); gsub(/NC_000016.9/, "16"); gsub(/NC_000017.10/, "17"); gsub(/NC_000018.9/, "18"); gsub(/NC_000019.9/, "19"); gsub(/NC_000020.10/, "20"); gsub(/NC_000021.8/, "21"); gsub(/NC_000022.10/, "22"); gsub(/NC_000023.10/, "23"); gsub(/NC_000024.9/, "24"); print;}' GCF_000001405.25.v4  > GCF_000001405.25.v5

Chromosome Y	CM000686.1	=	

awk '{print($1,$2,$3,$4,$5)}' GCF_000001405.25.v2 > GCF_000001405.25.v3

awk '{print $1":"$2":"$4":"$5"\t"$3}' < GCF_000001405.25.v5 > GCF_000001405.25.v6

grep -w -F -f cad_snps.txt GCF_000001405.25.v6 > GCF_000001405.25.v7

###############################################################
############## Match CAD variants with BIM files ##############
###############################################################

cad <- read.table("CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt", skip=15, header = TRUE)

files <- list.files(path="UKBB_G_Income/", pattern=".bim")

bim <- data.table()
for(i in files){
	temp <- fread(paste0("UKBB_G_Income/", i))
	temp$V2 <- paste0(temp$V1, ":", temp$V4, ":", temp$V5, ":", temp$V6)
	assign(gsub("\\..*","", i), temp)
# 	bim <- rbind(bim, temp)
}

files.2 <- gsub("\\..*","", files)
for(i in files.2){
	temp <- get(i)
	temp$V2 <- paste0(temp$V1, ":", temp$V4, ":", temp$V5, ":", temp$V6)
	assign(gsub("\\..*","", i), temp)
}

for(i in files.2){
	temp <- get(i)
	write.table(temp, paste0(i, ".bim"), sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

cad_ukbb <- inner_join(cad, bim, by = c("chr" = "V1", "position_hg19" = "V4", "A1"="V5", "A2" = "V6"))

# write.table(cad_ukbb[,c("V2", "effect_allele", "effect_weight")], "RM_CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt", sep ="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

######################################################
############## Update BIM file variants ##############
######################################################

files <- list.files(path="UKBB_G_Income/", pattern=".bim")

for(i in files){
	temp <- fread(paste0("UKBB_G_Income/", i))
	assign(gsub("\\..*","", i), temp)
}

files.2 <- gsub("\\..*","", files)

for(i in files.2){
	temp <- get(i)
	temp$V2 <- paste0(temp$V1, ":", temp$V4, ":", temp$V5, ":", temp$V6)
	assign(gsub("\\..*","", i), temp)
}

for(i in files.2){
	temp <- get(i)
	write.table(temp, paste0(i, ".bim"), sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

########################################################
################## Calculate Pred R^2 ##################
########################################################

library(data.table)
library(dplyr)
library(rcompanion)

# ukb <- read.table("gctb/MI/SAIGE/ukbb_prs_pheno_mi_test.txt", header = TRUE)
ukb <- read.table("gctb/MI/SAIGE/ukbb_prs_pheno_mi.txt", header=TRUE)

# files <- list.files(pattern=".profile")
# for(i in files){assign(gsub("\\..*","", i), fread(i))}

lky.chr1 <- read.table("cad_chr1.profile", header = T)
lky.all  <- lky.chr1
for (i in seq(2, 22)){
	print(paste("Doing", i))
	lky.chri <- read.table(paste0("cad_chr", i, ".profile"), header = T)
	lky.all$SCORE <- lky.all$SCORE + lky.chri$SCORE
}
lky.ukb <- inner_join(ukb, lky.all, by = c("ID" = "FID"))

resLocke <- data.frame(matrix(0, nrow = 1, ncol = 2))
resLocke[1, 1] <- "CAD full"
mod <- glm(lky.ukb$mi_flag ~ lky.ukb$SCORE, family="binomial")
resLocke[1, 2] <- nagelkerke(mod)$Pseudo.R.squared.for.model.vs.null[3]
print(resLocke)

colnames(resLocke) <- c("QC", "Pred_R2")
write.table(resLocke, "../cad_results.txt", col.names = T, row.names = F, quote = F, sep = " ")


#    --exclude UKBB_G_Income/LD_PLINK/LD_BED/bim${SLURM_ARRAY_TASK_ID}.dups 