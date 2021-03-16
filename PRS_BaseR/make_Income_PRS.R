##############################################################################
## Title: Income PRS
## Notes: Requires 48GB node
## Author: Regina Manansala
## Date Created: 05-March-2020
## Date Modified: 12-March-2020
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
dat <- fread("ukbb_prs_complete.txt", header=TRUE)

### Get SNP lists
snps <- fread("SNP_List/all_snps_w_chr.txt", header=TRUE)

start <- proc.time()
for(c in 1:22){
	### Subset dat for income SNPs. Take only half the sample
	h.snps <- subset(snps, Phenotype=="Income"  & Chromosome == c)
	if(nrow(h.snps) == 0){
		next
	}
	index1 <- which(colnames(dat) %in% as.character(h.snps$SNP_1))
	index2 <- which(colnames(dat) %in% as.character(h.snps$SNP_2))
	index <- union(index1, index2)
	print(length(index))
	h.dat <- dat[, c(1:21, index), with=FALSE]

	### No missing income values

	## Take a 50% sample from h.dat
	samp <- sample(1:315620, size=157810, replace=FALSE)
	samp2 <- setdiff(1:315620, samp)

	h.build <- h.dat[samp,]
	h.test <- h.dat[samp2,]

	### Clean up
# 	rm(list=c("dat", "snps", "h.dat"))
# 	gc()

	### Run a association tests for all these SNPs, save betas, in this 50% sample
	beta <- NULL
	p <- NULL
	y <- as.numeric(h.build$income)
	for(i in 1:length(index)){
		colnum <- 21+i
		x <- as.matrix(h.build[, ..colnum])[,1]
		mod <- lm(y ~ x)
		beta[i] <- summary(mod)$coef[2,1]
		p[i] <- summary(mod)$coef[2,4]
	}

	### Create PRS by summing beta*G (take 2-G when beta is negative, and (-1)*beta)

	## G matrix corresponding to positive betas
	if(length(which(beta>0)) == 0){
		G.positive <- NULL
	} else {
		G_index1 <- which(beta > 0)
		colnum <- 21 + G_index1
		G.positive <- as.matrix(h.test[, ..colnum])
	}
	
# 	G_index1 <- which(beta > 0)
# 	colnum <- 21 + G_index1
# 	G.positive <- as.matrix(h.test[, ..colnum])

	## G matrix corresponding to negative betas
	
	if(length(which(beta < 0)) == 0){
		G.negative <- NULL
		G.neg.v2 <- NULL
		beta.flip <- NULL
	} else {
		G_index2 <- which(beta < 0)
		colnum <- 21 + G_index2
		G.negative <- as.matrix(h.test[, ..colnum])
		G.neg.v2 <- 2 - G.negative
		beta.flip <- -1*beta[G_index2]
	}
	
# 	G_index2 <- which(beta < 0)
# 	colnum <- 21 + G_index2
# 	G.negative <- as.matrix(h.test[, ..colnum])
# 	G.neg.v2 <- 2 - G.negative
# 	beta.flip <- -1*beta[G_index2]

	## Cbind the G matrices
	G <- cbind(G.positive, G.neg.v2)

	## Clean up again
	rm(list=c("G.positive", "G.negative", "G.neg.v2", "h.build"))
	gc()
	
	if(exists('G_index1') == "FALSE" & is.null(beta.flip) == "FALSE"){
		beta.final <- beta.flip
	} 
	if(exists('G_index1') == "TRUE" & is.null(beta.flip) == "TRUE") {
		beta.final <- beta[G_index1]
	}
	if(exists('G_index1') == "TRUE" & is.null(beta.flip) == "FALSE"){
		beta.final <- c(beta[G_index1], beta.flip)
	}

	if(exists('G_index1') == "TRUE"){
		rm(G_index1)
		gc()
	}
	
	## Matrix multiplication of G*abs(beta) to get PRS
	prs <- G %*% beta.final

	## Histogram of PRS
	pdf(paste0("PRS/Income/income-prs-hist-chr", c, ".pdf"))
	hist(prs)
	dev.off()

	## Apply the PRS in the other 50% sample for income
	y <- as.numeric(h.test$income)
	x <- prs
	mod <- lm(y~x)
	print(summary(mod))

	## If it checks out, make a new file with just the PRS and phenotype and ID
	ID <- h.test$ID
	dat.out <- data.frame(ID, y, x)
	colnames(dat.out) <- c("ID", "income", "income_PRS")
# 	assign(paste0("dat_chr", c), dat.out)
	write.table(dat.out, paste0("PRS/Income/income-PRS-chr", c, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
}
end <- proc.time()

save.image(file = "PRS/Income/make_Income_prs.RData")





