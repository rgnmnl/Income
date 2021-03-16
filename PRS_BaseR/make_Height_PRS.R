##############################################################################
## Title: Height PRS
## Notes: Requires 48GB node
## Author: Regina Manansala
## Date Created: 10-March-2020
## Date Modified: 03-April-2020
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
dat <- fread("ukbb_prs_complete.txt", header=TRUE)

### Get SNP lists
snps <- fread("SNP_List/all_snps_w_chr.txt", header=TRUE)

###snps <- fread("SNP_List/all_snps_v2.txt", header=TRUE)

### Build PRS for Height by Chromosome

start <- proc.time()
for(c in 1:22){
	### Subset dat for Height SNPs. Take only half the sample
	h.snps <- subset(snps, Phenotype=="Height" & Chromosome == c)
	index1 <- which(colnames(dat) %in% as.character(h.snps$SNP_1))
	index2 <- which(colnames(dat) %in% as.character(h.snps$SNP_2))
	index <- union(index1, index2)
	print(length(index))
	h.dat <- dat[, c(1:21, index), with=FALSE]

	## Take a 50% sample from h.dat
	samp <- sample(1:315620, size=157810, replace=FALSE)
	samp2 <- setdiff(1:315620, samp)
	
	h.build <- h.dat[samp,]
	h.test <- h.dat[samp2,]

	### Clean up
# 	rm(list=c("dat", "snps", "h.dat"))
# 	gc()
	
	### Run a association tests for all these SNPs with height, save betas, in this 50% sample
	beta <- NULL
	p <- NULL
	y <- as.numeric(h.build$height)
	
	for(i in 1:length(index)){
		colnum <- 21+i
		x <- as.matrix(h.build[, ..colnum])[,1]
		mod <- lm(y ~ x)
		beta[i] <- summary(mod)$coef[2,1]
		p[i] <- summary(mod)$coef[2,4]
	}
	
	### Create PRS by summing beta*G (take 2-G when beta is negative, and (-1)*beta)

	## G matrix corresponding to positive betas
	G_index1 <- which(beta > 0)
	colnum <- 21 + G_index1
	G.positive <- as.matrix(h.test[, ..colnum])

	## G matrix corresponding to negative betas
	G_index2 <- which(beta < 0)
	colnum <- 21 + G_index2
	G.negative <- as.matrix(h.test[, ..colnum])
	G.neg.v2 <- 2 - G.negative
	beta.flip <- -1*beta[G_index2]

	## Cbind the G matrices
	G <- cbind(G.positive, G.neg.v2)
	assign(paste0("G_chr", c), G)
	
	## Clean up again
	rm(list=c("G.positive", "G.negative", "G.neg.v2", "h.build"))
	gc()

	beta.final <- c(beta[G_index1], beta.flip)

	## Matrix multiplication of G*abs(beta) to get PRS
	prs <- G %*% beta.final
	assign(paste0("prs_chr", c), prs)

	## Histogram of PRS
# 	pdf(paste0("PRS/Height/height-prs-hist-chr", c, ".pdf"))
# 	hist(prs)
# 	dev.off()
	
	## Apply the PRS in the other 50% sample for Height.
	y <- as.numeric(h.test$height)
	x <- prs
	mod <- lm(y~x)
	print(summary(mod))

	## If it checks out, make a new file with just the PRS and Height
	ID <- h.test$ID
	dat.out <- data.frame(ID, y, x)
	colnames(dat.out) <- c("ID", "height", "height_PRS")
	assign(paste0("dat_chr", c), dat.out)
# 	write.table(dat.out, paste0("PRS/Height/height-PRS-chr", c, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE)
}
end <- proc.time()

save.image(file = "PRS/Height/make_Height_prs.RData")

# for(c in 1:21){
# 
# chr_list <- h.snps %>% subset(., Chromosome == c) %>% select(., c("SNP_1", "SNP_2")) %>% gather(., "Variable", "SNP", c("SNP_1", "SNP_2")) %>% select(., "SNP") %>% as.list() %>% unlist()
# h.build.sub <- h.build %>% select(., c(1:21, which(names(h.build) %in% chr_list)))
# 
# beta <- NULL
# p <- NULL
# y <- as.numeric(h.build.sub$height)
# 
# for(i in 1:length(index)){
# colnum <- 21+i
# x <- as.matrix(h.build[, ..colnum])[,1]
# mod <- lm(y ~ x)
# beta[i] <- summary(mod)$coef[2,1]
# p[i] <- summary(mod)$coef[2,4]
# }
# }




