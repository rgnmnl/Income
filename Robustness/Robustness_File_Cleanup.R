##############################################################################
## Title: Robustness Reference File Data Cleanup
## Author: Regina Manansala
## Date Created: 14-August-2019
## Date Modified: 14-August-2019
##############################################################################

tyrrell <- fread("~/Documents/UKBB Income/BMI_Geno_Files/Tyrrell_2016_BMI_Pleio_SNPS.txt", sep = "", header = FALSE) %>% as.list() %>% unlist()
tyrrell_snp <- grep('rs', tyrrell, value = TRUE)
tyrrell_snp <- paste(unlist(t(tyrrell_snp)), collapse=" ")
tyrrell_snp <- gsub('BMI ', "", tyrrell_snp)
tyrrell_snp <- as.data.frame(tyrrell_snp)
tyrrell_snp <- separate_rows(tyrrell_snp, tyrrell_snp, sep = " ", convert = FALSE)
tyrrell_snp  <- grep("^rs*", tyrrell_snp$tyrrell_snp, value = TRUE)
write.table(tyrrell_snp[!(tyrrell_snp$tyrrell_snp %in% tyrrell_excl),], "Tyrrell_2016_BMI_Pleio_SNPS.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
tyrrell_excl<- c("rs11030104", "rs13107325", "rs3888190", "rs17001654", "rs2075650", "rs9925964", "rs2033529")


View(tyrrell_snp[!(tyrrell_snp$tyrrell_snp %in% tyrrell_excl),])

speliotes <- fread("~/Documents/UKBB Income/BMI_Geno_Files/Speliotes_2010_BMI_SNPs.txt", sep = "", header = FALSE) %>% as.list() %>% unlist()
spel_snp <- grep('rs', speliotes, value = TRUE)
spel_snp <- paste(unlist(t(spel_snp)), collapse=" ")
spel_snp <- as.data.frame(spel_snp)
spel_snp <- separate_rows(spel_snp, spel_snp, sep = " ", convert = FALSE)
spel_snp  <- grep("^rs*", spel_snp$spel_snp, value = TRUE)
write.table(spel_snp, "~/Documents/UKBB Income/BMI_Geno_Files/Speliotes_2010_BMI_SNPs.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

locke <- "rs657452 rs12286929 rs7903146 rs10132280 rs17094222 rs7599312 rs2365389 rs2820292 rs12885454 rs16851483 rs1167827 rs758747 rs1928295 rs9925964 rs11126666 rs2650492 rs6804842 rs4740619 rs13191362 rs3736485 rs17001654 rs11191560 rs1528435 rs1000940 rs2033529 rs11583200 rs9400239 rs10733682 rs11688816 rs11057405 rs11727676 rs3849570 rs6477694 rs7899106 rs2176598 rs2245368 rs17724992 rs7243357 rs2033732 rs9641123 rs7164727 rs492400 rs2080454 rs7239883 rs2836754 rs9914578 rs977747 rs9374842 rs4787491 rs1441264 rs17203016 rs16907751 rs13201877 rs9540493 rs1460676 rs6465468 rs6091540 rs7715256 rs2176040 rs1016287 rs4671328 rs758747 rs879620 rs12446632 rs11074446 rs6567160 rs17066842 rs9944545 rs11030104 rs10835210"
locke <- as.data.frame(locke)
locke_snp <- separate_rows(locke, locke, sep = " ", convert = FALSE)
