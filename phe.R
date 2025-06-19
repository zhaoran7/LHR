#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create phenotype data for GWAS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
dat0 <- readRDS(file = "/work/data/ukb/phe/ukb.phe.rds") # phenotype data extracted from UK Biobank

dat <- dat0 %>% filter(ethnic.c=="White") %>% rename(IID=eid) %>% mutate(FID = IID)
	dat$leg <- dat$height - dat$sitting_height
	dat$lhr = dat$leg / dat$height
	dat$lhr.res <- residuals(lm(lhr ~ height, data=dat, na.action=na.exclude)) # 

dat <- dat %>% dplyr::select(FID,IID,age,sex,height,leg,lhr,matches("^PC[0-9]+$"))
write.table(dat, "/work/data/ukb/phe/ukb.pheno", append=FALSE, quote=FALSE, row.names=FALSE)