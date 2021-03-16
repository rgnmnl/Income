##############################################################################
## Title: UKBB Depression GLM PRS
## Author: Regina Manansala
## Date Created: 13-August-2020
## Date Modified: 18-August-2020
##############################################################################

## List score files for test subset 
files <- list.files(pattern="*_test.txt")
files <- gsub("\\..*","", files)
files <- files[-1]

## Read test phenotype
dep_pheno_test <- fread("dep_pheno_test.txt")
dep_pheno_test <- dep_pheno_test[, 2:16]

## Read test subset scores
for(i in files){
	foo <- fread(paste0(i, ".txt")) %>% right_join(., dep_pheno_test, by = c("IID" = "ID"))
	assign(i, foo)
}

## Run regression of score vs. depression phenotype
prs_test_glm <- as.data.frame(matrix(0, ncol = 2, nrow = 22))
# prs_test_glm <- data.table(CHR=1:22, P=NA)
for(i in 1:length(files)){
# 	summary(glm(dep ~ SCORE1_AVG, data=get(files[i]), family = binomial(link="logit")))$coef[2,] %>% print()
	prs_test_glm[i, 1] <- files[i]
	pval <- summary(glm(dep ~ SCORE1_AVG, data=get(files[i]), family = binomial(link="logit")))$coef[2,4]
	prs_test_glm[i, 2] <- pval
}
prs_test_glm$V1 <- gsub("^dep_prs_|_test$", "", prs_test_glm$V1)
names(prs_test_glm) <- c("CHR", "P")

prs_test_glm$num <- gsub("chr", "", prs_test_glm$CHR) %>% as.numeric()
prs_test_glm <- prs_test_glm[order(prs_test_glm$num),]

## Export
write.table(prs_test_glm[, 1:2], "Depression_PRS_test_glm.txt", sep = "\t", row.names = FALSE, col.names=TRUE, quote =FALSE)

###################################################

## List PLINK score results
score <- list.files(pattern=".sscore")
score <- gsub("\\..*","", score)

## Read build phenotype
dep_pheno_build <- fread("dep_pheno_build.txt")
dep_pheno_build <- dep_pheno_build[, 2:16]

## Combine build and test pheno
dep_pheno <- rbind(dep_pheno_build, dep_pheno_test)

## Read score files
for(i in score){
	foo <- fread(paste0(i, ".sscore")) %>% right_join(., dep_pheno, by = c("IID" = "ID"))
	assign(i, foo)
}

## Run regression of score vs. depression phenotype
prs_glm <- as.data.frame(matrix(0, ncol = 2, nrow = 22))
for(i in 1:length(score)){
	summary(glm(dep ~ SCORE1_AVG, data=get(score[i]), family = binomial(link="logit")))$coef[2,] %>% print()
	prs_glm[i, 1] <- score[i]
	pval <- summary(glm(dep ~ SCORE1_AVG, data=get(score[i]), family = binomial(link="logit")))$coef[2,4]
	prs_glm[i, 2] <- pval
}
prs_glm$V1 <- gsub("^dep_prs_|_test$", "", prs_glm$V1)
names(prs_glm) <- c("CHR", "P")

prs_glm$num <- gsub("chr", "", prs_glm$CHR) %>% as.numeric()
prs_glm <- prs_glm[order(prs_glm$num),]

## Export
write.table(prs_glm[, 1:2], "Depression_PRS_glm.txt", sep = "\t", row.names = FALSE, col.names=TRUE, quote =FALSE)
