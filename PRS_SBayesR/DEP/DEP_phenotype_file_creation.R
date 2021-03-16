##############################################################################
## Title: UKBB Depression phenotype file creation
## Author: Regina Manansala
## Date Created: 07-August-2020
## Date Modified: 11-August-2020
##############################################################################

library(data.table)
library(dplyr)

## Read UKBB file
source("ukb43095.r")

## Get diagnosis variables and subset depression codes
ukb_dep <- select(bd, c("f.eid", "f.2010.0.0", "f.2090.0.0", "f.41202.0.0", "f.41202.0.1", "f.41204.0.0", "f.41204.0.1"))

ukb_dep$primary <- ifelse(grepl("F32|F33|F34|F38|F39", ukb_dep$f.41202.0.0), 1, 0)
ukb_dep$primary_v2 <- ifelse(grepl("F32|F33|F34|F38|F39", ukb_dep$f.41202.0.1), 1, 0)

ukb_dep$secondary <- ifelse(grepl("F32|F33|F34|F38|F39", ukb_dep$f.41204.0.0), 1, 0)
ukb_dep$secondary_v2 <- ifelse(grepl("F32|F33|F34|F38|F39", ukb_dep$f.41204.0.1), 1, 0)

ukb_dep$nerves <- ifelse(ukb_dep$f.2010.0.0 == "Yes", 1, 0)
ukb_dep$gp_dep <- ifelse(ukb_dep$f.2090.0.0 == "Yes", 1, 0)

ukb_dep$dep <- ukb_dep$primary + ukb_dep$primary_v2 + ukb_dep$secondary + ukb_dep$secondary_v2 + ukb_dep$nerves + ukb_dep$gp_dep
ukb_dep$dep <- ifelse(ukb_dep$dep >= 1, 1, 0)

## Read prs phenotype file and match with new depression data
ukbb_prs <- fread("../ukbb_prs_complete_v2.txt")
pheno <- ukbb_prs[,1:21]
geno <- ukbb_prs[,c(1, 22:8642)]

pheno <- left_join(pheno, ukb_dep[, c("f.eid", "dep")], by = c("ID" = "f.eid"))
pheno <- pheno[, c(1:20, 22)] %>% rename(., depression = dep.y)

ukbb_prs_v2 <- left_join(pheno, geno, by = "ID")

## Export
write.table(ukbb_prs_v2, "../ukbb_prs_complete_v2.txt", sep = "\t", row.names =FALSE, col.names=TRUE, quote = FALSE)

## Make build subset for GCTA
dep_build <- fread("dep_build_pheno_gcta_v2.txt")

write.table(pheno[pheno$ID %in% dep_build$V1, c(1,1,21)], "dep_build_pheno_gcta_v2.txt", sep = "\t", row.names =FALSE, col.names=FALSE, quote = FALSE)

############################################################

## Make build and test subsets for PRS
build <- fread("../PRS/PLINK_DEP_V2/dep_pheno_build.txt")
test <- fread("../PRS/PLINK_DEP_V2/dep_pheno_test.txt")

build_v2 <- pheno[pheno$ID %in% build$ID, c(1, 1:14, 21)]
names(build_v2) <- names(build)
test_v2 <- pheno[pheno$ID %in% test$ID, c(1, 1:14, 21)]
names(test_v2) <- names(test)

############################################################

## Get GCTA results
gcta[gcta$SNP %in% c("rs10127497", "rs6699744", "rs6424532", "rs7548151", "rs40465", "rs3132685", "rs112348907", "rs3807865", "rs2402273", "rs263575", "rs1021363", "rs10501696", "rs9530139", "rs28541419", "rs10929355", "rs5011432", "rs1554505"),]

## Get SNP info
pvar <- list.files(pattern="*.pvar", path="/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/")
# files <- list.files(pattern="^edunum.*\\.fastGWA$")
pvar

snps <- do.call(rbind, lapply(pvar, function(x) fread(paste0("/raid-05/SPH/pauer/UKBB/UKBB_G_Income/LD_PLINK/", x), stringsAsFactors = FALSE)))

snps[snps$ID %in% c("rs10127497", "rs6699744", "rs6424532", "rs7548151", "rs40465", "rs3132685", "rs112348907", "rs3807865", "rs2402273", "rs263575", "rs1021363", "rs10501696", "rs9530139", "rs28541419", "rs10929355", "rs5011432", "rs1554505"),]

############################################################

## Get plink regression results
files <- list.files(pattern="*.glm.logistic$")
# files <- list.files(pattern="^edunum.*\\.fastGWA$")
files

plink <- do.call(rbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE)))
pval_sub <- gcta %>% subset(., P < 5*10^-4) %>% arrange(., P)

plink[plink$ID %in% c("rs10127497", "rs6699744", "rs6424532", "rs7548151", "rs40465", "rs3132685", "rs112348907", "rs3807865", "rs2402273", "rs263575", "rs1021363", "rs10501696", "rs9530139", "rs28541419", "rs10929355", "rs5011432", "rs1554505"),]

############################################################

## Subset UKBB by depression medication 

library(data.table)
library(dplyr)

source('ukb23544_v2.r')
inc <- bd
source('ukb25323_v2.r')
gen <- bd
source('ukb29565_v2.r')
edu <- bd
source('ukb40060.r')
smk <- bd
source("ukb43095.r")
dep <- bd
source('ukb21755_v2.r')
# Self-reported Non-cancer illness: f.20002.0.0
	## Cases:  1289,  1291
	## Controls: 1289,  1291, 1286
# Treatment/medication: f.20003.0.0

ukbb <- left_join(dep, smk, by = "f.eid") %>% 
		left_join(., edu, by='f.eid') %>% 
		left_join(gen, by = 'f.eid') %>% 
		left_join(inc, by = 'f.eid') %>% 
		left_join(., bd[, c("f.21000.0.0", "f.31.0.0", "f.34.0.0", "f.52.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", "f.22009.0.1", "f.20003.0.0", "f.22009.0.2", "f.22011.0.0", "f.22012.0.0")], 
			by = c("f.21000.0.0", "f.31.0.0", "f.34.0.0", "f.52.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", "f.22009.0.1", "f.22009.0.2"))

ukbb$income[ukbb$f.738.0.0 == 'Less than 18,000'] <- 1
ukbb$income[ukbb$f.738.0.0 == '18,000 to 30,999'] <- 2
ukbb$income[ukbb$f.738.0.0 == '31,000 to 51,999'] <- 3
ukbb$income[ukbb$f.738.0.0 == '52,000 to 100,000'] <- 4
ukbb$income[ukbb$f.738.0.0 == 'Greater than 100,000'] <- 5

names(ukbb)[names(ukbb) == "f.eid"] <- "ID"
#names(ukbb)[names(ukbb) == "f.738.0.0"] <- "income"
names(ukbb)[names(ukbb) == "f.50.0.0"] <- "height" 
names(ukbb)[names(ukbb) == "f.31.0.0"] <- "sex" 
names(ukbb)[names(ukbb) == "f.22001.0.0"] <- "sex_genetic" 
names(ukbb)[names(ukbb) == "f.54.0.0"] <- "loc_code"
names(ukbb)[names(ukbb) == "f.34.0.0"] <- "birth_yr"
names(ukbb)[names(ukbb) == "f.52.0.0"] <- "birth_mo"
names(ukbb)[names(ukbb) == "f.6142.0.0"] <- "empl_status"
names(ukbb)[names(ukbb) == "f.21001.0.0"] <- "BMI"
names(ukbb)[names(ukbb) == "f.21003.0.0"] <- "Age"
names(ukbb)[names(ukbb) == "f.132.0.0"] <- "Job"
names(ukbb)[names(ukbb) == "f.22617.0.0"] <- "SOC_Code"
names(ukbb)[names(ukbb) == "f.1558.0.0"] <- "alc_int_freq"
names(ukbb)[names(ukbb) == "f.20160.0.0"] <- "ever_smoked"
names(ukbb)[names(ukbb) == "f.20414.0.0"] <- "alc_dr_freq"
names(ukbb)[names(ukbb) == "f.20544.0.1"] <- "mental_health"
names(ukbb)[names(ukbb) == "f.845.0.0"] <- "Age_Ed"

# names(ukbb)[42:51] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
names(ukbb)[300:309] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")


ukbb$diabetes <- ifelse(is.na(ukbb$f.2443.0.0) |ukbb$f.2443.0.0 == "No" | ukbb$f.2443.0.0 == "Prefer not to answer" | ukbb$f.2443.0.0 == "Do not know", 0, 1) #diabetes
ukbb$asthma <- ifelse(ukbb$f.3786.0.0 == -3 | ukbb$f.3786.0.0 == -1 | is.na(ukbb$f.3786.0.0), 0, 1)
ukbb$mi_flag <- ifelse(is.na(ukbb$f.42000.0.0), 0, 1) #myocardial infarction
ukbb$cancer_dx_dt <- ifelse(is.na(ukbb$f.40005.0.0), 0, 1) #cancer diagnosis date
ukbb$cancer_sr <- ifelse(is.na(ukbb$f.20001.0.0), 0, 1) #self-reported cancer

ukbb$race[ukbb$f.21000.0.0 == "British"] <- "British"
ukbb$race[ukbb$f.21000.0.0 == "Irish"] <- "Irish"
ukbb$race[ukbb$f.21000.0.0 == "Any other white background" | ukbb$f.21000.0.0 == "White"] <- "Other White"
ukbb$race[ukbb$f.21000.0.0 == "Chinese"] <- "Chinese"
ukbb$race[ukbb$f.21000.0.0 == "Caribbean" | ukbb$f.21000.0.0 == "African" | ukbb$f.21000.0.0 == "Any other Black background" | ukbb$f.21000.0.0 == "Black or Black British"] <- "Black"
ukbb$race[ukbb$f.21000.0.0 == "Mixed" | ukbb$f.21000.0.0 == "White and Black Caribbean" | ukbb$f.21000.0.0 == "White and Black African" | ukbb$f.21000.0.0 == "White and Asian" | ukbb$f.21000.0.0 == "Any other mixed background" | ukbb$f.21000.0.0 == "Other ethnic group" | ukbb$f.21000.0.0 == "Do not know" | ukbb$f.21000.0.0 == "Prefer not to answer" | ukbb$f.21000.0.0 == "Any other Asian background" | ukbb$f.21000.0.0 == "Asian or Asian British"] <- "Other"
ukbb$race[ukbb$f.21000.0.0 == "Indian" | ukbb$f.21000.0.0 == "Pakistani" | ukbb$f.21000.0.0 == "Bangladeshi"] <- "South Asian"

ukbb$location[ukbb$loc_code == 11012] <- "Barts"
ukbb$location[ukbb$loc_code == 11021] <- "Birmingham"
ukbb$location[ukbb$loc_code == 11011] <- "Bristol"
ukbb$location[ukbb$loc_code == 11008] <- "Bury"
ukbb$location[ukbb$loc_code == 11003] <- "Cardiff"
ukbb$location[ukbb$loc_code == 11024] <- "Cheadle"
ukbb$location[ukbb$loc_code == 11020] <- "Croydon"
ukbb$location[ukbb$loc_code == 11005] <- "Edinburgh"
ukbb$location[ukbb$loc_code == 11004] <- "Glasgow"
ukbb$location[ukbb$loc_code == 11018] <- "Hounslow"
ukbb$location[ukbb$loc_code == 11010] <- "Leeds"
ukbb$location[ukbb$loc_code == 11016] <- "Liverpool"
ukbb$location[ukbb$loc_code == 11001] <- "Manchester"
ukbb$location[ukbb$loc_code == 11017] <- "Middlesborough"
ukbb$location[ukbb$loc_code == 11009] <- "Newcastle"
ukbb$location[ukbb$loc_code == 11013] <- "Nottingham"
ukbb$location[ukbb$loc_code == 11002] <- "Oxford"
ukbb$location[ukbb$loc_code == 11007] <- "Reading"
ukbb$location[ukbb$loc_code == 11014] <- "Sheffield"
ukbb$location[ukbb$loc_code == 11006] <- "Stoke"
ukbb$location[ukbb$loc_code == 11022] <- "Swansea"
# ukbb$location[ukbb$loc_code == 11023] <- "Wrexham"
ukbb$location[ukbb$loc_code == 11025] <- "Cheadle (Imaging)"
ukbb$location[ukbb$loc_code == 11026] <- "Reading (Imaging)"
ukbb$location[ukbb$loc_code == 11027] <- "Newcastle (imaging)"

names(ukbb)[1714] <- "Relatedness_pairing"
names(ukbb)[1715] <- "Relatedness_factor"

ukbb <- subset(ukbb, is.na(Relatedness_factor))

# head(ukbb[which(grepl("F30|F31|F448|F2", ukbb$f.41202.0.0) | grepl("F30|F31|F448|F2", ukbb$f.41202.0.1) | grepl("F30|F31|F448|F2", ukbb$f.41204.0.0) | grepl("F30|F31|F448|F2", ukbb$f.41204.0.1)), c("f.41202.0.0", "f.41202.0.1", "f.41204.0.0", "f.41204.0.1")], 20)
# dim(ukbb[which(grepl("F30|F31|F448|F2", ukbb$f.41202.0.0) | grepl("F30|F31|F448|F2", ukbb$f.41202.0.1) | grepl("F30|F31|F448|F2", ukbb$f.41204.0.0) | grepl("F30|F31|F448|F2", ukbb$f.41204.0.1)), c("f.41202.0.0", "f.41202.0.1", "f.41204.0.0", "f.41204.0.1")])
# ukbb_noschiz <- subset(ukbb_rel, !(grepl("F30|F31|F448|F2", ukbb$f.41202.0.0) | grepl("F30|F31|F448|F2", ukbb$f.41202.0.1) | grepl("F30|F31|F448|F2", ukbb$f.41204.0.0) | grepl("F30|F31|F448|F2", ukbb$f.41204.0.1)))

ukbb_noschiz <- subset(ukbb, !(grepl("F30|F31|F448|F2", f.41202.0.0) | grepl("F30|F31|F448|F2", f.41202.0.1) | grepl("F30|F31|F448|F2", f.41204.0.0) | grepl("F30|F31|F448|F2", f.41204.0.1)))

ukbb_noillness <- subset(ukbb_noschiz, !(f.20002.0.0 %in% c("1291", "1289")))

rx_list <- c("1141202024", "1141153490", "1141195974", "1140867078", "1140867494", "1141171566", "2038459704", "1140872064", "1140879658", "1140867342", 
	"1140867420", "1140882320", "1140872216", "1140910358", "1141200458", "1141172838", "1140867306", "1140867180", "1140872200", "1140867210", "1140867398", 
	"1140882098", "1140867184", "1140867168", "1140863416", "1140909802", "1140867498", "1140867490", "1140910976", "1140867118", "1140867456", "1140928916", 
	"1140872268", "1140867134", "1140867208", "1140867218", "1140867572", "1140879674", "1140909804", "1140867504", "1140868170", "1140879746", "1141152848", 
	"1141177762", "1140867444", "1140867092", "1141152860", "1140872198", "1140867244", "1140868172", "1140867304", "1140872072", "1140879750", "1140868120", 
	"1140872214", "1141201792", "1140882100", "1141167976")
	
ukbb_norx <- subset(ukbb_noillness, !(f.20003.0.0 %in% rx_list))
dim(ukbb_noillness[ukbb_noillness$f.20003.0.0 %in% rx_list, ])

#### CASES

# dim(ukbb_norx[which(grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.0) | grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.1) | grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.0) | grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41204.0.1)), ])

ukbb_norx$primary <- ifelse(grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.0), 1, 0)
ukbb_norx$primary_v2 <- ifelse(grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.1), 1, 0)

ukbb_norx$secondary <- ifelse(grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41204.0.0), 1, 0)
ukbb_norx$secondary_v2 <- ifelse(grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41204.0.1), 1, 0)

# table(ukbb_norx$f.2010.0.0)
# table(ukbb_norx$f.2090.0.0)

ukbb_norx$nerves <- ifelse(ukbb_norx$f.2010.0.0 == "Yes", 1, 0)
ukbb_norx$gp <- ifelse(ukbb_norx$f.2090.0.0 == "Yes", 1, 0)

ukbb_norx$dep <- ukbb_norx$primary + ukbb_norx$primary_v2 + ukbb_norx$secondary + ukbb_norx$secondary_v2 + ukbb_norx$nerves + ukbb_norx$gp
ukbb_norx$dep <- ifelse(ukbb_norx$dep >= 1, 1, NA)

# ukbb_cases <- subset(ukbb_norx, !(grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.0) | grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.1) | grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41202.0.0) | grepl("F32|F33|F34|F38|F39", ukbb_norx$f.41204.0.1)) | f.2010.0.0 == "Yes" | f.2090.0.0 == "Yes")
ukbb_cases <- subset(ukbb_norx, dep == 1)

### CONTROLS

ukbb_norx$dep  <- ifelse(ukbb_norx$f.2010.0.0 == "No" & ukbb_norx$f.2090.0.0 == "No", 0, ukbb_norx$dep)

rx_list_controls <- c("1140867820", "1140867948", "1140879616", "1140867938", "1140867690", "1141190158", "1141151946", "1140921600", 
"1140879620", "1141201834", "1140867152", "1140909806", "1140879628", "1140867640", "1141200564", "1141151982", 
"1140916288", "1141180212", "1140867860", "1140867952", "1140879540", "1140867150", "1140909800", "1140867940", 
"1140879544", "1140879630", "1140867856", "1140867726", "1140867884", "1140867922", "1140910820", "1140879556", 
"1141152732", "1140867920", "1140882244", "1140867852", "1140867818", "1141174756", "1140867916", "1140867888", 
"1140867850", "1140867624", "1140867876", "1141151978", "1140882236", "1140867878", "1201", "1140882312", "1140867758", 
"1140867712", "1140867914", "1140867944", "1140879634", "1140867756", "1140867934", "1140867960", "1140916282", 
"1141200570", "1141152736")

ukbb_controls <- subset(ukbb_norx, dep == 0) %>% subset(., !(f.20003.0.0 %in% rx_list_controls)) %>%
	subset(., primary == 0 & primary_v2 == 0 & secondary == 0 & secondary_v2 == 0)

ukbb_comb <- rbind(ukbb_cases, ukbb_controls)

ukbb_white <- select(ukbb_comb, -starts_with("f.")) %>% subset(., race == "British")

	
# ukbb_rel <- left_join(ukbb, 
# 	bd[, c("f.eid", "f.21000.0.0", "f.31.0.0", "f.34.0.0", "f.52.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", "f.22009.0.1", "f.22009.0.2", "f.22011.0.0", "f.22012.0.0")], 
# 	by = c("race" = "f.21000.0.0", "sex" = "f.31.0.0", "birth_yr" = "f.34.0.0", "birth_mo" = "f.52.0.0", "f.53.0.0" = "f.53.0.0", "loc_code" = "f.54.0.0", "BMI" = "f.21001.0.0", "PC1" = "f.22009.0.1", "PC2" = "f.22009.0.2"))
