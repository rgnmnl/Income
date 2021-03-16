##############################################################################
## Title: Bootstrap
## Version: 1
## Author: Regina Manansala
## Date Created: 26-October-2020
## Date Modified: 06-November-2020
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
pheno.1 <-  as.character(args[1])
pheno.2 <- as.character(args[2])

library(data.table)
library(dplyr)
library(tibble)
library(stringr)

source("sBr_var_cov.R")

prop.var <- NULL
prop.cov <- NULL

test.ids <- fread("inc_prs_gctb/gctb/test_subset_IDs.txt") #%>% `[[`("V1") %>% as.list() %>% unlist()

# start <- proc.time()
start <- Sys.time()
for(i in 1:100000){
	sample.ids <- sample(test.ids$V1, size=length(test.ids$V1), replace=TRUE) %>% as.data.frame() %>% rename(., "V1" = `.`)
		
	prop.var[i] <- sbr_var(pheno.1, ids = TRUE)
	prop.cov[i] <- sbr_cov(pheno.1, pheno.2, ids = TRUE)
}
# end <- proc.time()
end <- Sys.time()

prop.var <- sort(prop.var)
prop.cov <- sort(prop.cov)

quantiles <- rbind(variance = quantile(prop.var, probs=c(0.05, 0.95)), covariance = quantile(prop.cov, probs=c(0.05, 0.95)))

write.table(prop.var, paste0("inc_prs_gctb/gctb/bootstrap/", pheno.1, "_var.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(prop.cov, paste0("inc_prs_gctb/gctb/bootstrap/", pheno.1, "_", pheno.2, "_cov.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(quantiles, paste0("inc_prs_gctb/gctb/bootstrap/", pheno.1, "_", pheno.2, "_quantiles.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

print(end-start)
print(quantiles)