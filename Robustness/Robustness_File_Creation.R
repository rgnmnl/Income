##############################################################################
## Title: Income Robustness File Creation
## Author: Regina Manansala
## Date Created: 22-July-2019
## Date Modified: 30-August-2019
##############################################################################

library(snpStats)
library(data.table)
library(dplyr)
library(BGData)

## Read SNPs from literature
tyrrell <- fread("Tyrrell_2016_BMI_SNPS.txt", header = F)
tyrrell_pleio <- fread("Tyrrell_2016_BMI_Pleio_SNPS.txt", header = F)
speliotes <- fread("Speliotes_2010_BMI_SNPs.txt", header = F)
locke <- fread("Locke_2015_BMI_SNPs.txt", header = F)
locke_pleio <- fread("Locke_Pleio_SNPs.txt", header = T)

##################
# geno_path <- "~/Data/UKBB_Genotypes" %>% print()
# 
# geno_fam <- list.files(path=geno_path, pattern='chr(\\d|\\d\\d)_PLINK_v2.pgen')
# geno_bim <- list.files(path=geno_path, pattern='chr(\\d|\\d\\d)_PLINK_v2.pvar')
# geno_bed <- list.files(path=geno_path, pattern='chr(\\d|\\d\\d)_PLINK_v2.psam')
# 
# ##make new number column
# 
# fam <- "ukb_filt.fam"
# bed <- "ukb_filt.bed"
# bim <- "ukb_filt.bim"
# 
# fam <- "ukb_filt_maf05.fam"
# bed <- "ukb_filt_maf05.bed"
# bim <- "ukb_filt_maf05.bim"
# 
# plink <- read.plink(bed, bim, fam)
# 
# geno <- as(t(plink$genotypes), "numeric") %>% as.matrix()
#################

## Read BED files
geno_path <- "~/Data/UKBB_PH/BED_BMI_Robustness/"
chr1 <- BEDMatrix(paste0(geno_path, "chr1_BMI_rbst.bed"))
chr2 <- BEDMatrix(paste0(geno_path, "chr2_BMI_rbst.bed"))
chr3 <- BEDMatrix(paste0(geno_path, "chr3_BMI_rbst.bed"))
chr4 <- BEDMatrix(paste0(geno_path, "chr4_BMI_rbst.bed"))
chr5 <- BEDMatrix(paste0(geno_path, "chr5_BMI_rbst.bed"))
chr6 <- BEDMatrix(paste0(geno_path, "chr6_BMI_rbst.bed"))
chr7 <- BEDMatrix(paste0(geno_path, "chr7_BMI_rbst.bed"))
chr8 <- BEDMatrix(paste0(geno_path, "chr8_BMI_rbst.bed"))
chr9 <- BEDMatrix(paste0(geno_path, "chr9_BMI_rbst.bed"))
chr10 <- BEDMatrix(paste0(geno_path, "chr10_BMI_rbst.bed"))
chr11 <- BEDMatrix(paste0(geno_path, "chr11_BMI_rbst.bed"))
chr12 <- BEDMatrix(paste0(geno_path, "chr12_BMI_rbst.bed"))
chr13 <- BEDMatrix(paste0(geno_path, "chr13_BMI_rbst.bed"))
chr14 <- BEDMatrix(paste0(geno_path, "chr14_BMI_rbst.bed"))
chr15 <- BEDMatrix(paste0(geno_path, "chr15_BMI_rbst.bed"))
chr16 <- BEDMatrix(paste0(geno_path, "chr16_BMI_rbst.bed"))
chr17 <- BEDMatrix(paste0(geno_path, "chr17_BMI_rbst.bed"))
chr18 <- BEDMatrix(paste0(geno_path, "chr18_BMI_rbst.bed"))
chr19 <- BEDMatrix(paste0(geno_path, "chr19_BMI_rbst.bed"))
chr20 <- BEDMatrix(paste0(geno_path, "chr20_BMI_rbst.bed"))
chr21 <- BEDMatrix(paste0(geno_path, "chr21_BMI_rbst.bed"))

wg <- ColumnLinkedMatrix(chr1, chr2, chr3, chr4, chr5, chr6, chr7,
						chr8, chr9, chr10, chr11, chr12, chr13, chr14,
						chr15, chr16, chr17, chr18, chr19, chr20, chr21)
geno <- as.matrix(wg) %>% as.data.frame()
geno$ID <- gsub("_.*","",rownames(geno)) %>% as.integer()
foo <- gsub("_.*", "", colnames(geno))
foo[43] <- paste0(foo[43], ".1")
names(geno) <- foo


#####################

## Working Directory: /raid-05/SPH/pauer/UKBB_PH/INC/SEM_FILES/TABLES

# fam <- "../GENO_FILES/sem_all_chrom_rn.fam"
# bed <- "../GENO_FILES/sem_all_chrom.bed"
# bim <- "../GENO_FILES/sem_all_chrom.bim"
# 
# plink <- read.plink(bed, bim, fam)
# 
# geno <- as(plink$genotypes, "numeric") %>% as.data.table()
# geno$ID <- plink$fam$member

bmi_disc <- fread("bmi_geno_analysis_w_ed.txt")
bmi_rep <- fread("bmi_geno_replicate_w_ed.txt")

bmi_comb <- rbind(bmi_disc, bmi_rep) %>% select(., c("ID", "BMI", "Age", "sex", "location", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "income", "Age_Ed"))

geno_tyrrell <- select(geno, c(which(colnames(geno) %in% tyrrell$V1), "ID"))
geno_tyrrell_pleio <- select(geno, c(which(colnames(geno) %in% tyrrell_pleio$V1), "ID"))
geno_speliotes <- select(geno, c(which(colnames(geno) %in% speliotes$V1), "ID"))
geno_locke <- select(geno, c(which(colnames(geno) %in% locke$V1), "ID"))
geno_locke_pleio <- select(geno, c(which(colnames(geno) %in% locke_pleio$Locke_Pleio_SNP), "ID"))

tyrrell_bmi <- left_join(bmi_comb, geno_tyrrell, by = "ID") #%>% subset(., !is.na(PC1))
tyrrell_pleio_bmi <- left_join(bmi_comb, geno_tyrrell_pleio, by = "ID") #%>% subset(., !is.na(PC1))
speliotes_bmi <- left_join(bmi_comb, geno_speliotes, by = "ID") #%>% subset(., !is.na(PC1))
locke_bmi <- left_join(bmi_comb, geno_locke, by = "ID") #%>% subset(., !is.na(PC1))
locke_pleio_bmi <- left_join(bmi_comb, geno_locke_pleio, by = "ID") #%>% subset(., !is.na(PC1))

write.table(tyrrell_bmi, "tyrrell_bmi.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(tyrrell_pleio_bmi, "tyrrell_pleio_bmi.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(speliotes_bmi, "speliotes_bmi.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(locke_bmi, "locke_bmi.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(locke_pleio_bmi, "locke_pleio_bmi.txt", sep="\t", row.names=F, col.names=T, quote=F)

######## IMPUTATION ##########

tyrrell_bmi_imp <- do.call(cbind, lapply(tyrrell_bmi, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x})) %>% as.data.table()
tyrrell_pleio_bmi_imp <- do.call(cbind, lapply(tyrrell_pleio_bmi, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x})) %>% as.data.table()
speliotes_bmi_imp <- do.call(cbind, lapply(speliotes_bmi, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x})) %>% as.data.table()
locke_bmi_imp <- do.call(cbind, lapply(locke_bmi, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x})) %>% as.data.table()
locke_pleio_bmi_imp <- do.call(cbind, lapply(locke_pleio_bmi, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x})) %>% as.data.table()

summary(tyrrell_pleio_bmi[, c(18:20)])
head(tyrrell_pleio_bmi_imp[, c(18:20)], 20)

summary(speliotes_bmi[, c(18:20)])
head(speliotes_bmi_imp[, c(18:20)], 20)

summary(locke_bmi[, c(18:20)])
head(locke_bmi_imp[, c(18:20)], 20)

summary(locke_pleio_bmi[, c(18:20)])
head(locke_pleio_bmi_imp[, c(18:20)], 20)

write.table(tyrrell_bmi_imp, "tyrrell_bmi_imp.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(tyrrell_pleio_bmi_imp, "tyrrell_pleio_bmi_imp.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(speliotes_bmi_imp, "speliotes_bmi_imp.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(locke_bmi_imp, "locke_bmi_imp.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(locke_pleio_bmi_imp, "locke_pleio_bmi_imp.txt", sep="\t", row.names=F, col.names=T, quote=F)

####### UNION AND INTERSECTION

locke_pleio <- locke_pleio[locke_pleio$Locke_Pleio_SNP != "", "Locke_Pleio_SNP"]
names(locke_pleio) <- "V1"
snps_union <- bind_rows(tyrrell, tyrrell_pleio, speliotes, locke, locke_pleio)
snps_union <- snps_union[unique(snps_union$V1),]

geno_union <- select(geno, c(which(colnames(geno) %in% snps_union$V1), "ID"))
union_bmi <- left_join(bmi_comb, geno_union, by = "ID") #%>% subset(., !is.na(PC1))
union_bmi_imp <- do.call(cbind, lapply(union_bmi, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x})) %>% as.data.table()

write.table(union_bmi_imp, "union_bmi_imp.txt", sep="\t", row.names=F, col.names=T, quote=F)


snps_intersect <- Reduce(intersect, list(as.list(tyrrell$V1),as.list(tyrrell_pleio$V1),as.list(speliotes$V1), as.list(locke$V1), as.list(locke_pleio$V1))) %>% unlist()

# tyrrell_lst <- as.list(tyrrell$V1) %>% unlist()
# tyrrell_pleio_lst <- as.list(tyrrell_pleio$V1) %>% unlist()
# speliotes_lst <- as.list(speliotes$V1) %>% unlist()
# locke_lst <- as.list(locke$V1) %>% unlist()
# locke_pleio_lst <- as.list(locke_pleio$V1) %>% unlist()
# snps_intersect <- Reduce(intersect, list(tyrrell_lst, tyrrell_pleio_lst))
# snps_intersect <- Reduce(intersect, list(tyrrell_lst, tyrrell_pleio_lst, speliotes_lst, locke_lst, locke_pleio_lst))

geno_int <- select(geno, c(which(colnames(geno) %in% snps_intersect), "ID"))
intersect_bmi <- left_join(bmi_comb, geno_int, by = "ID") #%>% subset(., !is.na(PC1))
intersect_bmi_imp <- do.call(cbind, lapply(intersect_bmi, function(x) {if(is.numeric(x)) ifelse(is.na(x),mean(x,na.rm=T),x) else x})) %>% as.data.table()

write.table(intersect_bmi_imp, "intersect_bmi_imp.txt", sep="\t", row.names=F, col.names=T, quote=F)



######## CHECK

tyrrell_og <- fread("tyrrell_bmi_imp.txt")
tyrrell_pleio_og <- fread("tyrrell_pleio_bmi_imp.txt")
speliotes_og <- fread("speliotes_bmi_imp.txt")
locke_og <- fread("locke_bmi_imp.txt")
locke_pleio_og <- fread("locke_pleio_bmi_imp.txt")

og_files <- ls(pattern="_og")

for(i in og_files){
print(i)
dim(i)
}




