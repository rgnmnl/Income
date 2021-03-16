##############################################################################
## Title: UKBB BMI Phenotype File Creation
## Author: Regina Manansala
## Date Created: 24-October-2018
## Date Modified: 24-October-2018
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(MASS)

## Import UKBB Income data
source('../ukb23544.r')
inc <- bd
source('../ukb25323.r')

## Import ID exclusion list
# ids <- as.list(read.table("w40458_20181016.txt", header=FALSE))
# bd[bd$f.eid %in% ids$V1, 1:5]
# ukbb_inc <- subset(bd, !(f.eid %in% ids$V1))

## Merge UKBB data and INC data
# foo <- ukbb_inc
ukbb_inc <- left_join(bd, inc, by='f.eid')

## Rename Variables (Income Data)
names(ukbb_inc)[names(ukbb_inc) == "f.eid"] <- "ID"
names(ukbb_inc)[names(ukbb_inc) == "f.738.0.0"] <- "income"
names(ukbb_inc)[names(ukbb_inc) == "f.50.0.0"] <- "height" 
names(ukbb_inc)[names(ukbb_inc) == "f.31.0.0"] <- "sex" 
names(ukbb_inc)[names(ukbb_inc) == "f.22001.0.0"] <- "sex_genetic" 
names(ukbb_inc)[names(ukbb_inc) == "f.54.0.0"] <- "loc_code"
names(ukbb_inc)[names(ukbb_inc) == "f.34.0.0"] <- "birth_yr"
names(ukbb_inc)[names(ukbb_inc) == "f.52.0.0"] <- "birth_mo"
names(ukbb_inc)[names(ukbb_inc) == "f.6142.0.0"] <- "empl_status"
names(ukbb_inc)[names(ukbb_inc) == "f.21001.0.0"] <- "BMI"
names(ukbb_inc)[names(ukbb_inc) == "f.21003.0.0"] <- "Age"
names(ukbb_inc)[names(ukbb_inc) == "f.132.0.0"] <- "Job"
names(ukbb_inc)[names(ukbb_inc) == "f.22617.0.0"] <- "SOC_Code"
names(ukbb_inc)[14:23] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

ukbb_inc$diabetes <- ifelse(is.na(ukbb_inc$f.2443.0.0) |ukbb_inc$f.2443.0.0 == "No" | ukbb_inc$f.2443.0.0 == "Prefer not to answer" | ukbb_inc$f.2443.0.0 == "Do not know", 0, 1) #diabetes
ukbb_inc$asthma <- ifelse(ukbb_inc$f.3786.0.0 == -3 | ukbb_inc$f.3786.0.0 == -1 | is.na(ukbb_inc$f.3786.0.0), 0, 1)
ukbb_inc$mi_flag <- ifelse(is.na(ukbb_inc$f.42000.0.0), 0, 1) #myocardial infarction
ukbb_inc$cancer_dx_dt <- ifelse(is.na(ukbb_inc$f.40005.0.0), 0, 1) #cancer diagnosis date
ukbb_inc$cancer_sr <- ifelse(is.na(ukbb_inc$f.20001.0.0), 0, 1) #self-reported cancer

ukbb_inc$race[ukbb_inc$f.21000.0.0 == "British"] <- "British"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Irish"] <- "Irish"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Any other white background" | ukbb_inc$f.21000.0.0 == "White"] <- "Other White"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Chinese"] <- "Chinese"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Caribbean" | ukbb_inc$f.21000.0.0 == "African" | ukbb_inc$f.21000.0.0 == "Any other Black background" | ukbb_inc$f.21000.0.0 == "Black or Black British"] <- "Black"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Mixed" | ukbb_inc$f.21000.0.0 == "White and Black Caribbean" | ukbb_inc$f.21000.0.0 == "White and Black African" | ukbb_inc$f.21000.0.0 == "White and Asian" | ukbb_inc$f.21000.0.0 == "Any other mixed background" | ukbb_inc$f.21000.0.0 == "Other ethnic group" | ukbb_inc$f.21000.0.0 == "Do not know" | ukbb_inc$f.21000.0.0 == "Prefer not to answer" | ukbb_inc$f.21000.0.0 == "Any other Asian background" | ukbb_inc$f.21000.0.0 == "Asian or Asian British"] <- "Other"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Indian" | ukbb_inc$f.21000.0.0 == "Pakistani" | ukbb_inc$f.21000.0.0 == "Bangladeshi"] <- "South Asian"

ukbb_inc_sub <- dplyr::select(ukbb_inc, -starts_with("f."))

## Remove NA income
inc_nomiss_inc <- subset(ukbb_inc_sub, !is.na(ukbb_inc_sub$income) & ukbb_inc_sub$income != "Prefer not to answer" & ukbb_inc_sub$income != "Do not know")
inc_nomiss_inc$income <- factor(inc_nomiss_inc$income)

## Remove NA employment
inc_nomiss_emp <- subset(inc_nomiss_inc, !(inc_nomiss_inc$empl_status %in% c("None of the above", "Prefer not to answer")))
inc_nomiss_emp$empl_status <- factor(inc_nomiss_emp$empl_status)

## Remove locations with <1000
inc_locex <- subset(inc_nomiss_emp, loc_code != 11023)

## Rename locations
inc_locex$location[inc_locex$loc_code == 11012] <- "Barts"
inc_locex$location[inc_locex$loc_code == 11021] <- "Birmingham"
inc_locex$location[inc_locex$loc_code == 11011] <- "Bristol"
inc_locex$location[inc_locex$loc_code == 11008] <- "Bury"
inc_locex$location[inc_locex$loc_code == 11003] <- "Cardiff"
inc_locex$location[inc_locex$loc_code == 11024] <- "Cheadle"
inc_locex$location[inc_locex$loc_code == 11020] <- "Croydon"
inc_locex$location[inc_locex$loc_code == 11005] <- "Edinburgh"
inc_locex$location[inc_locex$loc_code == 11004] <- "Glasgow"
inc_locex$location[inc_locex$loc_code == 11018] <- "Hounslow"
inc_locex$location[inc_locex$loc_code == 11010] <- "Leeds"
inc_locex$location[inc_locex$loc_code == 11016] <- "Liverpool"
inc_locex$location[inc_locex$loc_code == 11001] <- "Manchester"
inc_locex$location[inc_locex$loc_code == 11017] <- "Middlesborough"
inc_locex$location[inc_locex$loc_code == 11009] <- "Newcastle"
inc_locex$location[inc_locex$loc_code == 11013] <- "Nottingham"
inc_locex$location[inc_locex$loc_code == 11002] <- "Oxford"
inc_locex$location[inc_locex$loc_code == 11007] <- "Reading"
inc_locex$location[inc_locex$loc_code == 11014] <- "Sheffield"
inc_locex$location[inc_locex$loc_code == 11006] <- "Stoke"
inc_locex$location[inc_locex$loc_code == 11022] <- "Swansea"
# inc_locex$location[inc_locex$loc_code == 11023] <- "Wrexham"
inc_locex$location[inc_locex$loc_code == 11025] <- "Cheadle (Imaging)"
inc_locex$location[inc_locex$loc_code == 11026] <- "Reading (Imaging)"
inc_locex$location[inc_locex$loc_code == 11027] <- "Newcastle (imaging)"


## Remove NA race
inc_raceex <- subset(inc_locex, !is.na(race))
inc_whiteonly <- subset(inc_locex, race == "British")

## Limit Age to 25-65
inc_agesub <- subset(inc_whiteonly, Age >= 25 & Age <= 65)

## Create inc_analysis dataframe
inc_analysis <- inc_agesub

########

##Create continuous income variable
inc_analysis$inc_con[inc_analysis$income == 'Less than 18,000'] <- 1
inc_analysis$inc_con[inc_analysis$income == '18,000 to 30,999'] <- 2
inc_analysis$inc_con[inc_analysis$income == '31,000 to 51,999'] <- 3
inc_analysis$inc_con[inc_analysis$income == '52,000 to 100,000'] <- 4
inc_analysis$inc_con[inc_analysis$income == 'Greater than 100,000'] <- 5

################ OR ##################

## Create Binary Income Variable
# inc_analysis$inc_bin[inc_analysis$income == 'Less than 18,000' | inc_analysis$income == '18,000 to 30,999'] <- 1 # "Low"
# inc_analysis$inc_bin[inc_analysis$income == '52,000 to 100,000' | inc_analysis$income == 'Greater than 100,000'] <- 0 # "High"
# inc_analysis$inc_bin[inc_analysis$income == '31,000 to 51,999'] <- NA
# table(inc_analysis$inc_bin, exclude=NULL)
# inc_analysis$inc_bin <- as.numeric(inc_analysis$inc_bin)

########


inc_analysis$race_fct <- as.factor(inc_analysis$race)

inc_analysis$employment <- as.character(inc_analysis$empl_status)
inc_analysis$employment <- as.factor(inc_analysis$employment)

inc_analysis$sexfct <- as.character(inc_analysis$sex)
inc_analysis$sexfct <- as.factor(inc_analysis$sexfct)

inc_analysis$location <- as.factor(inc_analysis$location)

## Double-check variables to exclude NAs/unwanted categories
# table(inc_analysis$income, exclude=NULL)
table(inc_analysis$employment, exclude=NULL)
table(inc_analysis$sexfct, exclude=NULL)
table(inc_analysis$race_fct, exclude=NULL)
table(inc_analysis$location, exclude=NULL)

# table(inc_analysis$diabetes, exclude=NULL)
# table(inc_analysis$asthma, exclude=NULL)
# table(inc_analysis$mi_flag, exclude=NULL)
# table(inc_analysis$cancer_dx_dt, exclude=NULL)
# table(inc_analysis$cancer_sr, exclude=NULL)
# table(inc_analysis$birth_mo, exclude=NULL)


#######################

## Merge SOC Income Data
# inc_analysis$Job_SOC <- substr(inc_analysis$Job, 1, 4)
# 
# soc_income_2003 <- fread("SOC_gross_annual_pay.txt")
# soc_income_2003$Code <- as.numeric(soc_income_2003$Code)
# soc_income_2003 <- subset(soc_income_2003, !is.na(Code))
# 
# test <- left_join(inc_analysis, soc_income_2003, by=c("SOC_Code"="Code"))
# test <- test[, -c(19:20, 22, 24:34)]
# names(test)[19] <- paste0("soc_", names(test)[19])
# names(test)[20] <- paste0("soc_", names(test)[20])
# test2 <- left_join(test, soc_income_2003, by=c("Job_SOC"="Code"))
# test2 <- test2[, -c(21:22, 24, 26:36)]
# names(test2)[21] <- paste0("ukb_", names(test2)[21])
# names(test2)[22] <- paste0("ukb_", names(test2)[22])
# 
# dim(test2[is.na(test2$soc_Mean),])
# dim(test2[is.na(test2$ukb_Mean),])
# 
# test2$ukb_Mean <- gsub(",","",test2$ukb_Mean)
# test2$soc_Mean <- gsub(",","",test2$soc_Mean)
# 
# test2$ukb_Mean_num <- as.numeric(test2$ukb_Mean)
# test2$soc_Mean_num <- as.numeric(test2$soc_Mean)
# 
# test2$diff_Mean <- test2$ukb_Mean_num - test2$soc_Mean_num
# 
# head(test2[!is.na(test2$diff_Mean), c("income", "ukb_Mean_num", "soc_Mean_num", "diff_Mean")], 50)
# head(inc_analysis[!is.na(inc_analysis$SOC_Code) & inc_analysis$SOC_Code != inc_analysis$Job_SOC, ])
# 
# 
# test2$ukb_Mean_cat[test2$ukb_Mean_num < 18000] <- "Less than 18,000"
# test2$ukb_Mean_cat[test2$ukb_Mean_num >= 18000 & test2$ukb_Mean_num < 30999] <- "18,000 to 30,999"
# test2$ukb_Mean_cat[test2$ukb_Mean_num >= 31000 & test2$ukb_Mean_num < 51999] <- "31,000 to 51,999"
# test2$ukb_Mean_cat[test2$ukb_Mean_num >= 52000 & test2$ukb_Mean_num <= 100000] <- "52,000 to 100,000"
# test2$ukb_Mean_cat[test2$ukb_Mean_num > 100000] <- "Greater than 100,000"
# test2$ukb_Mean_cat[is.na(test2$ukb_Mean_num)] <- "NA"
# 
# dim(test2[test2$diff_Mean >= 18000 & test2$diff_Mean < 30999 & !is.na(test2$diff_Mean),])
# head(test2[test2$income == test2$ukb_Mean_cat, c("Job_SOC", "SOC_Code", "income", "diff_Mean_cat")], 50)
# 
# test2$diff_Mean_cat[test2$diff_Mean > 100000] <- "Greater than 100,000"
# dim(test2[test2$diff_Mean < 18000 & !is.na(test2$diff_Mean),])
# 
# summary(inc_analysis$SOC_Code)
# head(inc_analysis[!is.na(inc_analysis$SOC_Code), ])
# dim(inc_analysis[!is.na(inc_analysis$SOC_Code),])
# head(inc_analysis[!is.na(inc_analysis$SOC_Code),], 20)
# inc_analysis$Job_SOC <- substr(inc_analysis$Job, 1, 4)
# head(inc_analysis[!is.na(inc_analysis$SOC_Code),c("SOC_Code", "Job_SOC")], 20)
# head(inc_analysis[!is.na(inc_analysis$SOC_Code) & inc_analysis$SOC_Code != inc_analysis$Job_SOC, ])


###################


## Remove NA BMI
inc_bmiex <- subset(inc_analysis, !is.na(BMI))

########## FOR BINARY INCOME #################

## Remove NA binary Income
incbin_nomiss <- subset(inc_bmiex, !is.na(inc_bin))

## Keep desired covariates
bmi_cov <- incbin_nomiss[,c(1,1,31,22,34,30,4:13)]
names(bmi_cov)[1:6] <- c("#FID", "IID", "Income", "Age", "Sex", "Location")
bmi_cov$Income <- ifelse(bmi_cov$Income == 1, 2, 1)

## Random Sample 80% of data for analysis/20% for replication
# install.packages("ISLR")
library(ISLR)
smp_siz = floor(0.8*nrow(bmi_cov))
set.seed(123)
sample <- sample(seq_len(nrow(bmi_cov)),size = smp_siz)
bmi_analysis <- bmi_cov[sample,]
bmi_replicate <- bmi_cov[-sample,]

########## FOR CONTINUOUS INCOME #################

inc_cov <- inc_bmiex[,c(1,1,31,22,34,30,4:13)]
names(inc_cov)[1:6] <- c("#FID", "IID", "Income", "Age", "Sex", "Location")

#Random Sample 80% of data for analysis/20% for replication
# install.packages("ISLR")
library(ISLR)
smp_siz = floor(0.8*nrow(inc_cov))
set.seed(123)
sample <- sample(seq_len(nrow(inc_cov)),size = smp_siz)
inc_con_analysis <- inc_cov[sample,]
inc_con_replicate <- inc_cov[-sample,]

## BMI Sanity Check
bmi_ch16 <- incbin_nomiss[,c(1,1,21)]
names(bmi_ch16)[1:3] <- c("#FID", "IID", "BMI")

## Create BMI phenotype file

write.table(bmi_analysis, "bmi_analysis.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(bmi_ch16, "bmi_ch16.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(bmi_analysis[1:3], "bmi_pheno.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(bmi_replicate, "bmi_replicate.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(inc_con_analysis, "inc_con_analysis.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
write.table(inc_con_replicate, "inc_con_replicate.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
