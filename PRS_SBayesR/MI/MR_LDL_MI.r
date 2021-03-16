##############################################################################
## Title: UKBB MI LDL Mendelian Randomization analysis
## Author: Regina Manansala
## Date Created: 03-February-2021
## Date Modified: 10-February-2021
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)

## Read MI SAIGE results
mi_files <- list.files(path="MI_SAIGE", pattern=".SAIGE.bgen.txt")

mi <- data.table()
for(i in mi_files){
	temp <- fread(paste0("/raid-04/SPH/pauer/manansa2/inc_prs/MR_LDL_MI/MI_SAIGE/", i))
	mi <- rbind(mi, temp)
}

## Read MI GCTA results
# mi_files <- list.files(path="~/Data/manansa2/inc_prs/MR_LDL_MI/MI_TEST_GCTA/", pattern= ".fastGWA")
# 
# mi <- data.table()
# for(i in mi_files){
# 	temp <- fread(paste0("~/Data/manansa2/inc_prs/MR_LDL_MI/MI_TEST_GCTA/", i))
# 	mi <- rbind(mi, temp)
# }

## Read LDL GCTA results
ldl_files <- list.files(path="~/Data/manansa2/inc_prs/MR_LDL_MI/LDL_TEST_GCTA/", pattern= ".fastGWA")

ldl <- data.table()
for(i in ldl_files){
	temp <- fread(paste0("~/Data/manansa2/inc_prs/MR_LDL_MI/LDL_TEST_GCTA/", i))
	ldl <- rbind(ldl, temp)
}

## Get top 100 hits from MI (GCTA) and LDL results (from build set gwas)
mi_build <- fread("../../inc_prs_gctb/gctb/MI/GCTA/mi_gcta_top100.txt")
ldl_build <- fread("ldl_gcta_top100.txt")

## Get top 100 build SNPs from test set gwas results
mi_test <- mi %>% subset(., mi$rsid %in% mi_build$SNP) %>% filter(., !(rsid == "rs10738609" & Allele2 == "T"))
ldl_test <- ldl %>% subset(., SNP %in% ldl_build$SNP)

## Assign x, y and flip betas
x <- ldl_test$BETA
y <- mi_test$BETA

index <- which(x < 0)
x[index] <- -x[index]
y[index] <- -y[index]

## Run MR
MR <- lm(y~x)
MR_ivw <- lm(y~0+x)
  
print(summary(MR))
print(summary(MR_ivw))

## Plot results
p <- MR %>% ggplot(aes(x=x, y=y)) + geom_point() +
	geom_smooth(method="lm", formula = y~x,
                se = FALSE, size = 0.5,
                aes(color = "Egger")) +
    geom_smooth(method="lm", formula = y~0+x,
                se = FALSE, size = 0.5,
                aes(color = "IVW")) +
    theme_minimal() +
    ggtitle("MR LDL & MI") +
    ylab("MI BETA") + ##FLIP with X
    xlab("LDL BETA") + ##FLIP with Y
    scale_color_manual(name = "MR", 
                       labels = c(paste0("Egger", "\n", "(Slope=", format(MR$coefficients[2], digits = 2), ",\n", "P=", signif(summary(MR)$coef[2,4], 5), ")"), 
                                  paste0("IVW", "\n", "(Slope=", format(MR_ivw$coefficients[1], digits = 2), ",\n", "P=", signif(summary(MR_ivw)$coef[1,4], 5), ")")),
                       values = c("Egger" = "darkgreen",
                                  "IVW" = "darkorange"))

## Export plot
pdf("MR_LDL_MI_SAIGE.pdf")
print(p)
dev.off()          
           
## Plot Egger result only                       
p <- MR %>% ggplot(aes(x=x, y=y)) + geom_point() +
	geom_smooth(method="lm", formula = y~x, se = FALSE, size = 0.5, aes(color = "Egger")) +
    theme_minimal() +
    ggtitle("MR LDL & MI") +
    ylab("LDL BETA") +
    xlab("MI BETA") +
    scale_fill_identity(name = "MR", 
                       labels = c(paste0("Egger", "\n", "(Slope=", format(MR$coefficients[2], digits = 2), ",\n", "P=", signif(summary(MR)$coef[2,4], 5), ")"), 
                       values = c("Egger" = "darkgreen")))
  
## Export plot
pdf("MR_LDL_MI_Egger.pdf")
print(p)
dev.off()



