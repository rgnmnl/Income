##############################################################################
## Title: PRS Workflow Step 2
## Version: 1
## Purpose: Position filtering of fastGWA results using P-value
## Author: Regina Manansala
## Date Created: 18-May-2019
## Date Modified: 27-May-2020
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)

ht <- fread("height_gcta_results.txt")
pval_sub <- ht %>% subset(., P < 5*10^-4) %>% arrange(., P)

chr <- as.list(pval_sub$CHR) %>% unlist()
pos <- as.list(pval_sub$POS) %>% unlist()

pos_sub <- ht
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

save.image("Height_Pos_Excl_1.R")

for(i in 1:nrow(pos_sub)){
   pos <- pos_sub[i, "POS"]
   chr <- pos_sub[i, "CHR"]
   min_pos <- pos - 500000
   max_pos <- pos + 500000
   foo <- subset(pos_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, P != min(P) | duplicated(P))
   pos_sub <- pos_sub[!(pos_sub$id %in% drop_var$id),]
   if(i == nrow(pos_sub)) break
}

write.table(pos_sub[,1:14], "crp_pval_filt_2.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(pos_sub[,1:14], "il6_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)

write.table(pos_sub[,1:14], "tnfa_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)




// 
// drop_var <- data.frame()
// 
// for(i in 1:nrow(pval_sub)){
//    pos <- pval_sub[i, "POS"]
//    chr <- pval_sub[i, "CHR"]
//    min_pos <- pos - 250000
//    max_pos <- pos + 250000
//    pos_sub <- subset(pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
//    drop_var <- subset(pos_sub, p.value != min(p.value) | duplicated(p.value)) %>% rbind(., drop_var)
//    #pos_sub <- pos_sub[!(pos_sub$id %in% drop_var$id),]
//    if(i == nrow(pos_sub)) break
// }
