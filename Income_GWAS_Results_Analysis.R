##############################################################################
## Title: UKBB Income PLINK GWAS Results Analysis
## Version: 2
## Author: Regina Manansala
## Date Created: 01-February-2019
## Date Modified: 05-March-2019
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

### IMPORT TABLES
setwd("~/Documents/UKBB Income/Continuous_Income/")
bmi_snp <- fread("BMI_nature14177.csv")
t2d_snp <- fread("~/Documents/UKBB Income/MR Data/T2D_SNPs.csv")
cad_snp <- fread("~/Documents/UKBB Income/MR Data/CAD_SNPs.csv")

# bmi_snp_comp <- fread("BMI_nature14177.csv")
# foo <- left_join(bmi_snp, bmi_snp_comp, by="SNP")

inc_bmi <- fread("con_inc_bmi.txt")
inc_t2d <- fread("con_inc_t2d.txt")
inc_cad <- fread("con_inc_cad.txt")

names(inc_bmi)[1] <- "CHROM"
names(inc_t2d)[1] <- "CHROM"
names(inc_cad)[1] <- "CHROM"

### CREATE/RENAME ALLELE COLUMNS IN PUBLISHED DATA
cad_alleles <- str_split_fixed(cad_snp$Effect_Non_effect_allele, '/', 2)
cad_snp$REF <- cad_alleles[,1]
cad_snp$ALT1 <- cad_alleles[,2]

names(t2d_snp)[4] <- "REF"
names(t2d_snp)[5] <- "ALT1"

### JOIN PUB DATA & PLINK DATA
bmi <- left_join(inc_bmi, bmi_snp, by=c("ID"="SNP"))
t2d <- full_join(inc_t2d, t2d_snp, by=c("ID"="SNP"))
cad <- full_join(inc_cad, cad_snp, by=c("ID"="SNP"))

# CHECK NUMBER OF MISMATCHED ALLELES
dim(bmi[bmi$REF.x == bmi$REF.y & !is.na(bmi$REF.y),])
dim(t2d[t2d$REF.x == t2d$REF.y & !is.na(t2d$REF.y),])
dim(cad[cad$REF.x == cad$REF.y & !is.na(cad$REF.y),])

# CREATE TRANSFORMED BETA VARIABLE
bmi$BETA_2 <- ifelse(bmi$REF.x == bmi$REF.y & !is.na(bmi$REF.y), bmi$BETA.x, (-1)*bmi$BETA.x)

t2d$BETA_2 <- ifelse(t2d$REF.x == t2d$REF.y & !is.na(t2d$REF.y), t2d$BETA, (-1)*t2d$BETA)
t2d_no_out <- subset(t2d, b != 0.3059)

cad$logOR <- str_sub(cad$OR, 1, 4) %>% as.numeric() %>% log()
cad$BETA_2 <- ifelse(cad$REF.x == cad$REF.y & !is.na(cad$REF.y), cad$BETA, (-1)*cad$BETA)
cad_no_out <- subset(cad, logOR < 0.35065687)

head(bmi[bmi$REF.x != bmi$REF.y & !is.na(bmi$REF.y), c("REF.x", "ALT1.x", "BETA.x", "BETA_2", "REF.y", "ALT1.y")], 5)
head(t2d[t2d$REF.x != t2d$REF.y & !is.na(t2d$REF.y), c("REF.x", "ALT1.x", "BETA", "BETA_2", "REF.y", "ALT1.y")], 5)
head(cad[cad$REF.x != cad$REF.y & !is.na(cad$REF.y), c("REF.x", "ALT1.x", "BETA", "BETA_2", "REF.y", "ALT1.y")], 5)

dim(bmi[bmi$BETA.x == bmi$BETA_2 & !is.na(bmi$BETA_2),])
dim(t2d[t2d$BETA == t2d$BETA_2 & !is.na(t2d$BETA_2),])
dim(cad[cad$BETA == cad$BETA_2 & !is.na(cad$BETA_2),])

## HEIGHT
plink_ht <- fread("Height/plink_ht.txt")
plink_fem_ht <- fread("Height/plink_fem_ht.txt")
plink_male_ht <- fread("Height/plink_male_ht.txt")
#pub_ht <- fread("Height/Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt")
# height <- left_join(plink_ht, pub_ht, by=c("ID"="SNP"))
# height$BETA_2 <- ifelse(height$REF == height$Tested_Allele & !is.na(height$Tested_Allele), height$BETA.x, (-1)*height$BETA.x)

#height <- subset(height, P.y < 5^(-10))

## PLINK2 BMI

# inc_bmi_plink2 <- fread("inc_bmi_plink2.txt")
# bmi_plink2 <- left_join(inc_bmi, inc_bmi_plink2, by="ID")
# dim(bmi_plink2[bmi_plink2$REF.x != bmi_plink2$REF.y & !is.na(bmi_plink2$REF.y),])
plink_bmi <- fread("BMI/plink_bmi.txt")
plink_fem_bmi <- fread("BMI/plink_female_bmi.txt")
plink_male_bmi <- fread("BMI/plink_male_bmi.txt")

## PLINK2 MI
plink_mi <- fread("../MI/plink_mi.txt")
plink_fem_mi <- fread("../MI/plink_female_mi.txt")
plink_male_mi <- fread("../MI/plink_male_mi.txt")


### MR ANALYSIS FUNCTION

MR_fxn <- function(df, y_beta, x_beta, dis){
  y <- df[[y_beta]]
  x <- df[[x_beta]]
  index <- which(x < 0)
  x[index] <- -x[index]
  y[index] <- -y[index]
  
  MR <- lm(y~x)
  MR_ivw <- lm(y~0+x)
  
  print(summary(MR))
  print(summary(MR_ivw))
  
  p <- MR %>% ggplot(aes(x=x, y=y)) + geom_point() +
    geom_smooth(method="lm", formula = y~x,
                se = FALSE, size = 0.5,
                aes(color = "Egger")) +
    geom_smooth(method="lm", formula = y~0+x,
                se = FALSE, size = 0.5,
                aes(color = "IVW")) +
    theme_minimal() +
    ggtitle(paste("MR", dis, "SNPs", sep=" ")) +
    ylab("Income Betas") + ##FLIP with X
    xlab(paste(dis, "Betas", sep=" ")) + ##FLIP with Y
    scale_color_manual(name = "MR", 
                       labels = c(paste0("Egger", "\n", "(Slope=", format(MR$coefficients[2], digits = 2), ",\n", "P=", signif(summary(MR)$coef[2,4], 5), ")"), 
                                  paste0("IVW", "\n", "(Slope=", format(MR_ivw$coefficients[1], digits = 2), ",\n", "P=", signif(summary(MR_ivw)$coef[1,4], 5), ")")),
                       values = c("Egger" = "darkgreen",
                                  "IVW" = "darkorange"))
  
  pdf(paste0("plot", gsub(" ", "", dis, fixed = TRUE), ".pdf"))
  print(p)
  dev.off()
}

# MR_fxn(df=bmi, y_beta="BETA_2", x_beta="BETA.y", dis="BMI")
# 
# MR_fxn(df=t2d, y_beta = "BETA_2", x_beta = "b", dis = "Type 2 Diabetes")
# MR_fxn(df=t2d_no_out, y_beta = "BETA_2", x_beta = "b", dis = "Type 2 Diabetes")
# 
# MR_fxn(df=cad, y_beta = "BETA_2", x_beta = "logOR", dis = "Coronary Artery Disease")
# MR_fxn(df=cad_no_out, y_beta = "BETA_2", x_beta = "logOR", dis = "Coronary Artery Disease")
# 
# MR_fxn(df=bmi_plink2, y_beta = "BETA.x", x_beta = "BETA.y", dis = "BMI Plink2")

MR_fxn(df=plink_ht, y_beta = "inc_BETA", x_beta = "ht_BETA", dis = "Height Plink2")
MR_fxn(df=plink_fem_ht, y_beta = "inc_BETA", x_beta = "ht_BETA", dis = "Height Female Plink2")
MR_fxn(df=plink_male_ht, y_beta = "inc_BETA", x_beta = "ht_BETA", dis = "Height Male Plink2")

MR_fxn(df=plink_bmi, y_beta = "inc_BETA", x_beta = "bmi_BETA", dis = "BMI Plink2")
MR_fxn(df=plink_fem_bmi, y_beta = "inc_BETA", x_beta = "bmi_BETA", dis = "BMI Female Plink2")
MR_fxn(df=plink_male_bmi, y_beta = "inc_BETA", x_beta = "bmi_BETA", dis = "BMI Male Plink2")

MR_fxn(df=plink_mi, y_beta = "inc_BETA", x_beta = "mi_logOR", dis = "MI Plink2")
MR_fxn(df=plink_fem_mi, y_beta = "inc_BETA", x_beta = "mi_logOR", dis = "MI Female Plink2")
MR_fxn(df=plink_male_mi, y_beta = "inc_BETA", x_beta = "mi_logOR", dis = "MI Male Plink2")

MR_fxn(df=plink_mi[plink_mi$mi_ID != "rs3798220",], y_beta = "inc_BETA", x_beta = "mi_logOR", dis = "MI Plink2 No Outlier")
MR_fxn(df=plink_fem_mi[plink_fem_mi$mi_ID != "rs3798220",], y_beta = "inc_BETA", x_beta = "mi_logOR", dis = "MI Female Plink2 No Outlier")
MR_fxn(df=plink_male_mi[plink_male_mi$mi_ID != "rs3798220",], y_beta = "inc_BETA", x_beta = "mi_logOR", dis = "MI Male Plink2 No Outlier")




