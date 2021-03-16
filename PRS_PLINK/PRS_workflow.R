##############################################################################
## Title: PRS Workflow Step 1
## Purpose: Create phenotype and covar files for fastGWA analysis
## Author: Regina Manansala
## Date Created: 07-August-2019
## Date Modified: 29-May-2020
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

source("/raid-05/SPH/pauer/UKBB_PH/ukb42080.r")

pheno <- fread("/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/GWAS/build_pheno.txt")
pheno_test <- fread("/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/GWAS/test_pheno.txt")

alc_build_pheno <- pheno[, c(1:2, 16)]
write.table(alc_build_pheno, "alc_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

smk_build_pheno <- pheno[, c(1:2, 17)]
write.table(smk_build_pheno, "smk_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

ht_build_pheno <- pheno[, c(1:2, 19)]
write.table(ht_build_pheno, "ht_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

dep_build_pheno <- pheno[, c(1:2, 18)]
write.table(dep_build_pheno, "dep_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

bmi_build_pheno <- pheno[, c(1:2, 20)]
write.table(bmi_build_pheno, "bmi_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

inc_build_pheno <- pheno[, c(1:2, 21)]
write.table(inc_build_pheno, "inc_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

inc_build_pheno <- pheno[, c(1:2, 21)]
write.table(inc_build_pheno, "inc_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

pheno_edu <- left_join(pheno, bd[, 1:2], by = c("IID" = "f.eid")) %>% 
	rename(., Edu = f.6138.0.0) %>%
	mutate(., Edu_bin = case_when(
			Edu %in% c("None of the above", "Prefer not to answer") ~ NA_real_,
			Edu %in% "College or University degree" ~ 1,
			Edu %in% c("A levels/AS levels or equivalent", "O levels/GCSEs or equivalent", "CSEs or equivalent", "NVQ or HND or HNC or equivalent", "Other professional qualifications eg: nursing, teaching") ~ 0,
			TRUE ~ as.numeric(Edu)),
		Edu_num = case_when(
			Edu ==  "College or University degree" ~ 20,
			Edu ==  "NVQ or HND or HNC or equivalent" ~ 19,
			Edu ==  "Other professional qualifications eg: nursing, teaching" ~ 15,
			Edu ==  "A levels/AS levels or equivalent" ~ 13,
			Edu %in%  c("O levels/GCSEs or equivalent", "CSEs or equivalent") ~ 10,
			Edu == "None of the above" ~ 7,
			Edu == "Prefer not to answer" ~ NA_real_,
			TRUE ~ as.numeric(Edu))
			)
	
pheno_test_edu <- inner_join(pheno_test, bd[, 1:2], by = c("IID" = "f.eid")) %>%
	rename(., Edu = f.6138.0.0) %>%
	mutate(., Edu_bin = case_when(
			Edu %in% c("None of the above", "Prefer not to answer") ~ NA_real_,
			Edu %in% "College or University degree" ~ 1,
			Edu %in% c("A levels/AS levels or equivalent", "O levels/GCSEs or equivalent", "CSEs or equivalent", "NVQ or HND or HNC or equivalent", "Other professional qualifications eg: nursing, teaching") ~ 0,
			TRUE ~ as.numeric(Edu)),
		Edu_num = case_when(
			Edu ==  "College or University degree" ~ 20,
			Edu ==  "NVQ or HND or HNC or equivalent" ~ 19,
			Edu ==  "Other professional qualifications eg: nursing, teaching" ~ 15,
			Edu ==  "A levels/AS levels or equivalent" ~ 13,
			Edu %in%  c("O levels/GCSEs or equivalent", "CSEs or equivalent") ~ 10,
			Edu == "None of the above" ~ 7,
			Edu == "Prefer not to answer" ~ NA_real_,
			TRUE ~ as.numeric(Edu))
		)

write.table(pheno_edu[!is.na(pheno_edu$Edu_bin), c(1,1,24)], "edubin_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_edu[!is.na(pheno_edu$Edu_num), c(1,1,25)], "edunum_build_pheno_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

pheno_edu_nomiss <- fread("edubin_build_pheno_gcta.txt")
covar_edubin <- subset(covar, V1 %in% pheno_edu_nomiss$V1)
qcovar_edubin <- subset(qcovar, V1 %in% pheno_edu_nomiss$V1)
write.table(covar_edubin, "edubin_build_covar_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(qcovar_edubin, "edubin_build_qcovar_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

pheno_edu_nomiss <- fread("edunum_build_pheno_gcta.txt")
covar_edunum <- subset(covar, V1 %in% pheno_edu_nomiss$V1)
qcovar_edunum <- subset(qcovar, V1 %in% pheno_edu_nomiss$V1)
write.table(covar_edunum, "edunum_build_covar_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(qcovar_edunum, "edunum_build_qcovar_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


write.table(pheno_edu[, c(1:15, 19)], "PLINK_Height/ht_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:15, 19)], "PLINK_Height/ht_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(pheno_edu[, c(1:15, 20)], "PLINK_BMI/bmi_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:15, 20)], "PLINK_BMI/bmi_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(pheno_edu[, c(1:15, 24)], "PLINK_EDU_bin/edubin_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:15, 24)], "PLINK_EDU_bin/edubin_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(pheno_edu[, c(1:15, 25)], "PLINK_EDU_num/edunum_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:15, 25)], "PLINK_EDU_num/edunum_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(pheno_edu[, c(1:15, 21)], "PLINK_Income/inc_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:15, 21)], "PLINK_Income/inc_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(pheno_edu[, c(1:16)], "PLINK_ALC/alc_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:16)], "PLINK_ALC/alc_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(pheno_edu[, c(1:15, 17)], "PLINK_DEP/dep_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:15, 17)], "PLINK_DEP/dep_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(pheno_edu[, c(1:15, 18)], "PLINK_SMK/smk_pheno_build.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(pheno_test_edu[, c(1:15, 18)], "PLINK_SMK/smk_pheno_test.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


################################################################################################################

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)
library(tibble)
library(stringr)
library(qqman)

files <- list.files(pattern="*.fastGWA")
# files <- list.files(pattern="^edunum.*\\.fastGWA$")
files

gcta <- do.call(rbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE)))
pval_sub <- gcta %>% subset(., P < 5*10^-4) %>% arrange(., P)

chr <- as.list(pval_sub$CHR) %>% unlist()
pos <- as.list(pval_sub$POS) %>% unlist()

pos_sub <- pval_sub[order(pval_sub$CHR, pval_sub$POS, pval_sub$P),]

for(i in 1:nrow(pos_sub)){
   pos <- pos_sub[i, "POS"]
   chr <- pos_sub[i, "CHR"]
   min_pos <- pos - 500000
   max_pos <- pos + 500000
   foo <- subset(pos_sub, CHR == chr & POS %between% c(min_pos,max_pos))
#    print(nrow(foo))
   drop_var <- subset(foo, P != min(P) | duplicated(P))
   #print(dim(drop_var))
   pos_sub <- pos_sub[!(pos_sub$SNP %in% drop_var$SNP),]
   if(i == nrow(pos_sub)) break
}

write.table(pos_sub, "DEP_Index_SNPS_v2.txt", sep = "\t", row.names=FALSE, col.names = TRUE, quote = FALSE)

table1 <- gcta[gcta$SNP %in% c("rs10127497", "rs6699744", "rs6424532", "rs7548151", "rs40465", "rs3132685", "rs112348907", "rs3807865", "rs2402273", "rs263575", "rs1021363", "rs10501696", "rs9530139", "rs28541419", "rs10929355", "rs5011432", "rs1554505"),]
write.table(table1, "broad_depression_table1_SNPs.txt", sep = "\t", row.names=FALSE, col.names = TRUE, quote = FALSE)

jpeg("dep_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(gcta, chr="CHR", bp="POS", p="P", cex = 0.6)
dev.off()

jpeg("dep_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(gcta$P)
dev.off()

