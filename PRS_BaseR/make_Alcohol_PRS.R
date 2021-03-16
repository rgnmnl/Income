##############################################################################
## Title: Alcohol PRS
## Notes: Requires 48GB node
## Author: Regina Manansala
## Date Created: 30-January-2020
## Date Modified: 30-January-2020
##############################################################################

library(data.table)
dat <- fread("ukbb_prs_complete.txt", header=TRUE)

### Get SNP lists
#snps <- fread("SNP_List/SNP_reference.txt", header=TRUE)

snps <- fread("SNP_List/all_snps_v2.txt", header=TRUE)

### Subset dat for SNPs. Take only half the sample
h.snps <- subset(snps, Phenotype=="Alcohol")
index1 <- which(colnames(dat) %in% as.character(h.snps$SNP_1))
index2 <- which(colnames(dat) %in% as.character(h.snps$SNP_2))
index <- union(index1, index2)
h.dat <- dat[, c(1:21, index), with=FALSE]

## Take a 50% sample from h.dat
samp <- sample(1:315620, size=157810, replace=FALSE)
samp2 <- setdiff(1:315620, samp)

h.build <- h.dat[samp,]
h.test <- h.dat[samp2,]

### Clean up
rm(list=c("dat", "snps", "h.dat"))
gc()


### Run a association tests for all these SNPs, save betas, in this 50% sample
beta <- NULL
p <- NULL
y <- as.numeric(h.build$alc)
for(i in 1:length(index)){
colnum <- 21+i
x <- as.matrix(h.build[, ..colnum])[,1]
if(length(unique(x))==1){
	beta[i] <- 0
	p[i] <- 1}else{

	mod <- lm(y ~ x)
	beta[i] <- summary(mod)$coef[2,1]
	p[i] <- summary(mod)$coef[2,4]}
}

### Create PRS by summing beta*G (take 2-G when beta is negative, and (-1)*beta)

## G matrix corresponding to positive betas
index1 <- which(beta >= 0)
colnum <- 21 + index1
G.positive <- as.matrix(h.test[, ..colnum])

## G matrix corresponding to negative betas
index2 <- which(beta < 0)
colnum <- 21 + index2
G.negative <- as.matrix(h.test[, ..colnum])
G.neg.v2 <- 2 - G.negative
beta.flip <- -1*beta[index2]

## Cbind the G matrices
G <- cbind(G.positive, G.neg.v2)

## Clean up again
rm(list=c("G.positive", "G.negative", "G.neg.v2", "h.build"))
gc()

beta.final <- c(beta[index1], beta.flip)

## Matrix multiplication of G*abs(beta) to get PRS
prs <- G %*% beta.final

## Histogram of PRS
pdf("~/Data/Mike-Patrick-PRS/alcohol-prs-hist.pdf")
hist(prs)
dev.off()

## Apply the PRS in the other 50% sample.
y <- as.numeric(h.test$alc)
x <- prs
mod <- lm(y~x)

## If it checks out, make a new file with just the PRS and phenotype
ID <- h.test$ID
dat.out <- data.frame(ID, y, x)
colnames(dat.out) <- c("ID", "alcohol", "alc_PRS")
write.table(dat.out, "~/Data/Mike-Patrick-PRS/alcohol-PRS.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)









