##############################################################################
## Title: UKBB Myocardial Infarction Analysis
## Author: Regina Manansala
## Date Created: 28-February-2019
## Date Modified: 05-March-2019
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(MASS)

#Import UKBB Income data
source('../ukb23544.r')
inc <- bd
source('../ukb25323.r')

#Merge UKBB data and INC data
ukbb_mi <- left_join(bd, inc, by='f.eid')

#Rename Variables (Income Data)
names(ukbb_mi)[names(ukbb_mi) == "f.eid"] <- "ID"
names(ukbb_mi)[names(ukbb_mi) == "f.738.0.0"] <- "income"
names(ukbb_mi)[names(ukbb_mi) == "f.50.0.0"] <- "height" 
names(ukbb_mi)[names(ukbb_mi) == "f.31.0.0"] <- "sex" 
names(ukbb_mi)[names(ukbb_mi) == "f.22001.0.0"] <- "sex_genetic" 
names(ukbb_mi)[names(ukbb_mi) == "f.54.0.0"] <- "loc_code"
names(ukbb_mi)[names(ukbb_mi) == "f.34.0.0"] <- "birth_yr"
names(ukbb_mi)[names(ukbb_mi) == "f.52.0.0"] <- "birth_mo"
names(ukbb_mi)[names(ukbb_mi) == "f.6142.0.0"] <- "empl_status"
names(ukbb_mi)[names(ukbb_mi) == "f.21001.0.0"] <- "BMI"
names(ukbb_mi)[names(ukbb_mi) == "f.21003.0.0"] <- "Age"
names(ukbb_mi)[names(ukbb_mi) == "f.132.0.0"] <- "Job"
names(ukbb_mi)[names(ukbb_mi) == "f.22617.0.0"] <- "SOC_Code"
names(ukbb_mi)[14:23] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

#Create binary disease variables
ukbb_mi$diabetes <- ifelse(is.na(ukbb_mi$f.2443.0.0) |ukbb_mi$f.2443.0.0 == "No" | ukbb_mi$f.2443.0.0 == "Prefer not to answer" | ukbb_mi$f.2443.0.0 == "Do not know", 0, 1) #diabetes
ukbb_mi$asthma <- ifelse(ukbb_mi$f.3786.0.0 == -3 | ukbb_mi$f.3786.0.0 == -1 | is.na(ukbb_mi$f.3786.0.0), 0, 1)
ukbb_mi$mi_flag <- ifelse(is.na(ukbb_mi$f.42000.0.0), 0, 1) #myocardial infarction
ukbb_mi$cancer_dx_dt <- ifelse(is.na(ukbb_mi$f.40005.0.0), 0, 1) #cancer diagnosis date
ukbb_mi$cancer_sr <- ifelse(is.na(ukbb_mi$f.20001.0.0), 0, 1) #self-reported cancer

#Create race variable
ukbb_mi$race[ukbb_mi$f.21000.0.0 == "British"] <- "British"
ukbb_mi$race[ukbb_mi$f.21000.0.0 == "Irish"] <- "Irish"
ukbb_mi$race[ukbb_mi$f.21000.0.0 == "Any other white background" | ukbb_mi$f.21000.0.0 == "White"] <- "Other White"
ukbb_mi$race[ukbb_mi$f.21000.0.0 == "Chinese"] <- "Chinese"
ukbb_mi$race[ukbb_mi$f.21000.0.0 == "Caribbean" | ukbb_mi$f.21000.0.0 == "African" | ukbb_mi$f.21000.0.0 == "Any other Black background" | ukbb_mi$f.21000.0.0 == "Black or Black British"] <- "Black"
ukbb_mi$race[ukbb_mi$f.21000.0.0 == "Mixed" | ukbb_mi$f.21000.0.0 == "White and Black Caribbean" | ukbb_mi$f.21000.0.0 == "White and Black African" | ukbb_mi$f.21000.0.0 == "White and Asian" | ukbb_mi$f.21000.0.0 == "Any other mixed background" | ukbb_mi$f.21000.0.0 == "Other ethnic group" | ukbb_mi$f.21000.0.0 == "Do not know" | ukbb_mi$f.21000.0.0 == "Prefer not to answer" | ukbb_mi$f.21000.0.0 == "Any other Asian background" | ukbb_mi$f.21000.0.0 == "Asian or Asian British"] <- "Other"
ukbb_mi$race[ukbb_mi$f.21000.0.0 == "Indian" | ukbb_mi$f.21000.0.0 == "Pakistani" | ukbb_mi$f.21000.0.0 == "Bangladeshi"] <- "South Asian"

#Create continuous income variable
ukbb_mi$inc_con[ukbb_mi$income == 'Less than 18,000'] <- 1
ukbb_mi$inc_con[ukbb_mi$income == '18,000 to 30,999'] <- 2
ukbb_mi$inc_con[ukbb_mi$income == '31,000 to 51,999'] <- 3
ukbb_mi$inc_con[ukbb_mi$income == '52,000 to 100,000'] <- 4
ukbb_mi$inc_con[ukbb_mi$income == 'Greater than 100,000'] <- 5

#Flag Non-Cancer Illnesses
ukbb_mi$nci_mi <- rowSums(ukbb_mi[,178:264] == 1075, na.rm=TRUE)
table(ukbb_mi$nci_mi)
ukbb_mi$nci_mi <- ifelse(ukbb_mi$nci_mi > 0, 1, 0)

ukbb_mi$nci_is <- rowSums(ukbb_mi[,178:264] == 1583, na.rm=TRUE)
table(ukbb_mi$nci_is)

ukbb_mi$nci_str <- rowSums(ukbb_mi[,178:264] == 1081, na.rm=TRUE)
table(ukbb_mi$nci_str)
ukbb_mi$nci_str <- ifelse(ukbb_mi$nci_str > 0, 1, 0)

ukbb_mi$nci_sh <- rowSums(ukbb_mi[,178:264] == 1086, na.rm=TRUE)
table(ukbb_mi$nci_sh)
ukbb_mi$nci_sh <- ifelse(ukbb_mi$nci_sh > 0, 1, 0)

ukbb_mi$nci_bh <- rowSums(ukbb_mi[,178:264] == 1491, na.rm=TRUE)
table(ukbb_mi$nci_bh)
ukbb_mi$nci_bh <- ifelse(ukbb_mi$nci_bh > 0, 1, 0)

######### THIS INFORMATION IS NOT AVAILABLE IN UPDATED UKBB DATA #########
#Count how many hypertensive patients are in data from self-reported Non-Cancer Illness
# hypertension <- c(1065, 1072, 1073)
# ukbb_inc$hyp_events <- rowSums(ukbb_inc[,178:264] == hypertension, na.rm=TRUE)
# ukbb_inc$hyp_flag <- ifelse(ukbb_inc$hyp_events > 0, 1, 0)
# table(ukbb_inc$hyp_flag)


######### THIS INFORMATION IS NOT RELEVANT TO MI ANALYSIS #########
# Remove NA income
# inc_nomiss_inc <- subset(ukbb_inc_sub, !is.na(ukbb_inc_sub$income) & ukbb_inc_sub$income != "Prefer not to answer" & ukbb_inc_sub$income != "Do not know")
# inc_nomiss_inc$income <- factor(inc_nomiss_inc$income)

# Remove NA employment
# inc_nomiss_emp <- subset(inc_nomiss_inc, !(inc_nomiss_inc$empl_status %in% c("None of the above", "Prefer not to answer")))
# inc_nomiss_emp$empl_status <- factor(inc_nomiss_emp$empl_status)

#Remove locations with <1000
# inc_locex <- subset(inc_nomiss_emp, loc_code != 11023)

#Rename locations
# inc_locex$location[inc_locex$loc_code == 11012] <- "Barts"
# inc_locex$location[inc_locex$loc_code == 11021] <- "Birmingham"
# inc_locex$location[inc_locex$loc_code == 11011] <- "Bristol"
# inc_locex$location[inc_locex$loc_code == 11008] <- "Bury"
# inc_locex$location[inc_locex$loc_code == 11003] <- "Cardiff"
# inc_locex$location[inc_locex$loc_code == 11024] <- "Cheadle"
# inc_locex$location[inc_locex$loc_code == 11020] <- "Croydon"
# inc_locex$location[inc_locex$loc_code == 11005] <- "Edinburgh"
# inc_locex$location[inc_locex$loc_code == 11004] <- "Glasgow"
# inc_locex$location[inc_locex$loc_code == 11018] <- "Hounslow"
# inc_locex$location[inc_locex$loc_code == 11010] <- "Leeds"
# inc_locex$location[inc_locex$loc_code == 11016] <- "Liverpool"
# inc_locex$location[inc_locex$loc_code == 11001] <- "Manchester"
# inc_locex$location[inc_locex$loc_code == 11017] <- "Middlesborough"
# inc_locex$location[inc_locex$loc_code == 11009] <- "Newcastle"
# inc_locex$location[inc_locex$loc_code == 11013] <- "Nottingham"
# inc_locex$location[inc_locex$loc_code == 11002] <- "Oxford"
# inc_locex$location[inc_locex$loc_code == 11007] <- "Reading"
# inc_locex$location[inc_locex$loc_code == 11014] <- "Sheffield"
# inc_locex$location[inc_locex$loc_code == 11006] <- "Stoke"
# inc_locex$location[inc_locex$loc_code == 11022] <- "Swansea"
# # inc_locex$location[inc_locex$loc_code == 11023] <- "Wrexham"
# inc_locex$location[inc_locex$loc_code == 11025] <- "Cheadle (Imaging)"
# inc_locex$location[inc_locex$loc_code == 11026] <- "Reading (Imaging)"
# inc_locex$location[inc_locex$loc_code == 11027] <- "Newcastle (imaging)"


#Remove NA race
mi_whiteonly <- subset(ukbb_mi, !is.na(race)) %>% subset(., race == "British")

#Limit Age to 25-65
mi_agesub <- subset(mi_whiteonly, Age >= 25 & Age <= 65)

#Remove NA BMI
# inc_bmiex <- subset(inc_analysis, !is.na(BMI))

##Remove duplicate relatedness factor - GENETIC RELATEDNESS NOT AVAILABLE IN INCOME DATA
##https://stackoverflow.com/questions/24011246/deleting-rows-that-are-duplicated-in-one-column-based-on-the-condition\
##s-of-anoth?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
# no.rel <- inc_agesub[order(inc_agesub[,'f.22011.0.0']),]
# no.rel.dup <- no.rel[!duplicated(no.rel$f.22011.0.0,incomparables=NA),]
# table(no.rel.dup$mi_flag)

##Remove individuals with NA PC1-10
mi_pccomp <- subset(mi_agesub, !is.na(mi_agesub$PC1))
table(mi_pccomp$mi_flag)

#Create subset with only CASES
mi_cases <- filter(mi_pccomp, mi_flag == 1)
dim(mi_cases)

#Create subset with only CONTROLS - NO MI AND NO STROKE
mi_controls <- filter(mi_pccomp, mi_flag == 0)
dim(mi_controls)
table(mi_controls$mi_flag)

##Remove individuals with self-reported stroke/ischaemic stroke
mi_nonci <- filter(mi_controls, nci_str == 0 & nci_is == 0 & nci_sh == 0 & nci_bh == 0)
table(mi_nonci$nci_str)
table(mi_nonci$nci_is)
table(mi_nonci$nci_sh)
table(mi_nonci$nci_bh)

## FOR MI
##Total Controls = 53500; Total Males = 41730; Total Females = 11770
#Create Male only, no stroke  dataset
mi_male <- mi_nonci[mi_nonci$sex == "Male",]
set.seed(111)
mi_male.samp <-  mi_male[sample(nrow(mi_male), size=41730, replace=FALSE),]
#Create Female only, no stroke dataset
mi_female <- mi_nonci[mi_nonci$sex == "Female",]
set.seed(222)
mi_fem.samp <-  mi_female[sample(nrow(mi_female), size=11770, replace=FALSE),]
#Combine Sampled Controls and Cases
mi_analysis <- bind_rows(mi_fem.samp, mi_male.samp, mi_cases)

########
#Convert variables to factors or characters

mi_analysis$race_fct <- as.factor(mi_analysis$race)
mi_analysis$employment <- as.character(mi_analysis$empl_status)
mi_analysis$employment <- as.factor(mi_analysis$employment)
mi_analysis$sexfct <- as.character(mi_analysis$sex)
mi_analysis$sexfct <- as.factor(mi_analysis$sexfct)
#mi.analysis$location <- as.factor(mi.analysis$location)

#Double-check variables to exclude NAs/unwanted categories

# table(inc_analysis$income, exclude=NULL)
table(mi_analysis$employment, exclude=NULL)
table(mi_analysis$sexfct, exclude=NULL)
table(mi_analysis$race_fct, exclude=NULL)
#table(inc_analysis$location, exclude=NULL)


#Create Covariate File
#Drop unnecessary variables
mi_cov <- dplyr::select(mi_analysis, -starts_with("f."))

#Checks
table(mi_cov$mi_flag)

#Covariate Files - IS
mi_cov1 <- mi_cov[, c(1,1,26,22,38,4:13)]
names(mi_cov1)[1:5] <- c("#FID", "IID", "MI", "Age", "Sex")
mi_cov1$MI <- ifelse(mi_cov1$MI == 1, 2, 1)


write.table(mi_cov1, "mi_cov.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(mi_cov1[mi_cov1$Sex == "Female", -5],  "mi_fem_cov.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(mi_cov1[mi_cov1$Sex == "Male", -5],  "mi_male_cov.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

## IMPORT RESULTS INTO R

fmi_dir <- "~/Data/UKBB_PH/INC/MI_ANALYSIS/MI_FEMALE/"
mmi_dir <- "~/Data/UKBB_PH/INC/MI_ANALYSIS/MI_MALE/"

mi_regression_files <- list.files(pattern="mi_chr(\\d|\\d\\d)_assoc.MI.glm.logistic$")
mi_results <- do.call(rbind, lapply(mi_regression_files, function(x) fread(x, stringsAsFactors = FALSE)))
names(mi_results)[1] <- "CHROM"
colnames(mi_results) <- paste("mi", colnames(mi_results), sep="_")
mi_snps <- read.csv("~/Data/UKBB_PH/INC/CAD_SNPs.csv", header=T)
plink_mi <- left_join(mi_results[mi_ID %in% mi_snps$SNP, c(1:5, 12, 15)], inc_results[inc_ID %in% mi_snps$SNP, c("inc_ID", "inc_BETA", "inc_P")], by=c("mi_ID"="inc_ID"))
plink_mi$mi_logOR <- log(plink_mi$mi_OR)


mi_regression_files <- list.files(path=fmi_dir, pattern="mi_chr(\\d|\\d\\d)_assoc.MI.glm.logistic$")
mi_female <- do.call(rbind, lapply(mi_regression_files, function(x) fread(paste0(fmi_dir, x), stringsAsFactors = FALSE)))
names(mi_female)[1] <- "CHROM"
colnames(mi_female) <- paste("mi", colnames(mi_female), sep="_")
mi_snps <- read.csv("~/Data/UKBB_PH/INC/CAD_SNPs.csv", header=T)
plink_mi_female <- left_join(mi_female[mi_ID %in% mi_snps$SNP, c(1:5, 12, 15)], inc_results[inc_ID %in% mi_snps$SNP, c("inc_ID", "inc_BETA", "inc_P")], by=c("mi_ID"="inc_ID"))
plink_mi_female$mi_logOR <- log(plink_mi_female$mi_OR)


mi_regression_files <- list.files(path=mmi_dir, pattern="mi_chr(\\d|\\d\\d)_assoc.MI.glm.logistic$")
mi_male <- do.call(rbind, lapply(mi_regression_files, function(x) fread(paste0(mmi_dir, x), stringsAsFactors = FALSE)))
names(mi_male)[1] <- "CHROM"
colnames(mi_male) <- paste("mi", colnames(mi_male), sep="_")
mi_snps <- read.csv("~/Data/UKBB_PH/INC/CAD_SNPs.csv", header=T)
plink_mi_male <- left_join(mi_male[mi_ID %in% mi_snps$SNP, c(1:5, 12, 15)], inc_results[inc_ID %in% mi_snps$SNP, c("inc_ID", "inc_BETA", "inc_P")], by=c("mi_ID"="inc_ID"))
plink_mi_male$mi_logOR <- log(plink_mi_male$mi_OR)


write.table(plink_mi, "plink_mi.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(plink_mi_female, "plink_female_mi.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(plink_mi_male, "plink_male_mi.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)






