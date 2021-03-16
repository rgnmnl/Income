##############################################################################
## Title: PRS Workflow Step 2
## Version: 2
## Purpose: Position filtering of fastGWA results using P-value
## Author: Regina Manansala
## Date Created: 27-May-2019
## Date Modified: 1-June-2020
##############################################################################


library(data.table)
library(dplyr)
library(tidyr)

inc <- fread("inc_gcta_results.txt")
pval_sub <- inc %>% subset(., P < 5*10^-4) %>% arrange(., P)

chr <- as.list(pval_sub$CHR) %>% unlist()
pos <- as.list(pval_sub$POS) %>% unlist()

pos_sub <- inc
n_excl <- 0
start <- proc.time()
for(i in 1:24660){
	min_pos <- pos[i] - 500000
	max_pos <- pos[i] + 500000
	n_excl <- n_excl + nrow(pos_sub[pos_sub$POS > min_pos & pos_sub$POS < max_pos & pos_sub$CHR == chr[i],])
	print(nrow(pos_sub[pos_sub$POS > min_pos & pos_sub$POS < max_pos & pos_sub$CHR == chr[i],]))
	pos_sub <- subset(pos_sub, !(CHR == chr[i] & POS > min_pos & POS < max_pos))
}
end <- proc.time()

save.image("Income_Pos_Excl_1.R")


pos_sub_2 <- pval_sub[order(pval_sub$CHR, pval_sub$POS, pval_sub$P),]

for(i in 1:nrow(pos_sub_2)){
   pos <- pos_sub_2[i, "POS"]
   chr <- pos_sub_2[i, "CHR"]
   min_pos <- pos - 500000
   max_pos <- pos + 500000
   foo <- subset(pos_sub_2, CHR == chr & POS %between% c(min_pos,max_pos))
#    print(nrow(foo))
   drop_var <- subset(foo, P != min(P) | duplicated(P))
   #print(dim(drop_var))
   pos_sub_2 <- pos_sub_2[!(pos_sub_2$SNP %in% drop_var$SNP),]
   if(i == nrow(pos_sub_2)) break
}
write.table(pos_sub_2, "BMI_Index_SNPS.txt", sep = "\t", row.names=FALSE, col.names = TRUE, quote = FALSE)


# library(data.table)
# library(dplyr)
# library(tidyr)
# 
# bmi <- fread("bmi_gcta_results.txt")
# pval_sub <- bmi %>% subset(., P < 5*10^-4) %>% arrange(., P)
# 
# chr <- as.list(pval_sub$CHR) %>% unlist()
# pos <- as.list(pval_sub$POS) %>% unlist()
# 
# pos_sub <- bmi
# n_excl <- 0
# start <- proc.time()
# for(i in 1:83715){
# 	min_pos <- pos[i] - 500000
# 	max_pos <- pos[i] + 500000
# 	n_excl <- n_excl + nrow(pos_sub[pos_sub$POS > min_pos & pos_sub$POS < max_pos & pos_sub$CHR == chr[i],])
# 	print(nrow(pos_sub[pos_sub$POS > min_pos & pos_sub$POS < max_pos & pos_sub$CHR == chr[i],]))
# 	pos_sub <- subset(pos_sub, !(CHR == chr[i] & POS > min_pos & POS < max_pos))
# }
# end <- proc.time()
# 
# save.image("BMI_Pos_Excl_1.R")

