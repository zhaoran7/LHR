#---------------------------------------
# Mendelian randomization (GSMR method)
#---------------------------------------
install.packages(c('survey'))
install.packages(devtools::install_github("jianyanglab/gsmr2"))
pacman::p_load(data.table,stringi,TwoSampleMR,plink,gsmr2,dplyr)

# LHR and T2D association as example
# Load X 
dat.X.raw <- fread("/data/GWAS/height/LHR.gz");dat.X.raw$P <- as.numeric(dat.X.raw$P)
  dat.X.sig <- dat.X.raw %>% dplyr::filter(P<=5e-08)
  dat.X.sig <- dat.X.sig %>% dplyr::filter(!(CHR == 6 & POS >= 28477797 & POS <= 33448354))
  dat.X.sig <- dat.X.sig %>% mutate(F_stat = (BETA^2) / (SE^2)) %>% dplyr::filter(F_stat >= 10)
  dat.X.top <- fread('/data/GWAS/height/lhr.jma.cojo')
  dat.X <- data.frame(dat.X.sig[SNP %in% dat.X.top$SNP, ])%>% 
    format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P',samplesize_col = 'N', eaf_col = 'EAF') 
# Load Y (FinnGen)
dat.Y.raw <- fread("/data/GWAS/CMD/T2D.gz")  
dat.Y <- data.frame(dat.Y.raw %>% merge(dat.X, by="SNP")) %>% 
  format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P', samplesize_col = 'N',eaf_col = 'EAF') 

# harmonised data
dat.CMD <- harmonise_data(dat.X, dat.Y, action=1) 

# GSMR format
gsmr_CMD <- dat.CMD %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, samplesize.exposure,
         beta.outcome, eaf.outcome, se.outcome, pval.outcome, se.exposure, pval.exposure, samplesize.outcome)
colnames(gsmr_CMD) <- c("SNP", "a1", "a2", "bzx", "bzx_n","bzy", "a1_freq", "bzy_se","bzy_pval", "bzx_se", "bzx_pval","bzy_n")
#save genetic variants and ea
write_tsv(gsmr_CMD, "/data/gsmr/gsmr_CMD.tsv")
write.table(gsmr_CMD[,c(1,2)],"/data/gsmr/gsmr_CMD_snps.allele", col.names=F, row.names=F, quote=F)
#linux command use GCTA software
#gcta64 --bfile /data/1kG/EUR --extract /data/gsmr/gsmr_CMD_snps.allele --update-ref-allele /data/gsmr/gsmr_CMD_snps.allele --recode --out /data/gsmr/gsmr_CMD
#LD matrix
snp_coeff_id = scan("/data/gsmr/gsmr_CMD.xmat.gz", what="", nlines=1)
snp_coeff = read.table("/data/gsmr/gsmr_CMD.xmat.gz", header=F, skip=2)
snp_id = Reduce(intersect, list(gsmr_CMD$SNP, snp_coeff_id))
gsmr_data = gsmr_CMD[match(snp_id, gsmr_CMD$SNP),]
snp_order = match(snp_id, snp_coeff_id)
snp_coeff_id = snp_coeff_id[snp_order]
snp_coeff = snp_coeff[, snp_order]
ldrho = cor(snp_coeff)
colnames(ldrho) = rownames(ldrho) = snp_coeff_id
dim(ldrho)
#gsmr data
gsmr_CMD <- read.table("/data/gsmr/gsmr_CMD.tsv",header=T) %>% filter(SNP %in% snp_id)
gsmr_CMD_xmat <- read.table("/data/gsmr/gsmr_CMD.xmat.gz",header = T)
#coeff col
snpfreq = gsmr_CMD$a1_freq             # allele frequencies of the SNPs
bzx = gsmr_CMD$bzx     # effects of the instruments on risk factor
bzx_se = gsmr_CMD$bzx_se       # standard errors of bzx
bzx_n = gsmr_CMD$bzx_n          # GWAS sample size for the risk factor
std_zx = std_effect(snpfreq, bzx, bzx_se, bzx_n)    # perform standardisation
gsmr_CMD$std_bzx = std_zx$b    # standardized bzx
gsmr_CMD$std_bzx_se = std_zx$se    # standardized bzx_se
head(gsmr_CMD)
#gsmr
bzx = gsmr_CMD$std_bzx
bzx_se = gsmr_CMD$std_bzx_se
bzx_pval = gsmr_CMD$bzx_pval
bzy = gsmr_CMD$bzy
bzy_se = gsmr_CMD$bzy_se
bzy_pval = gsmr_CMD$bzy_pval
n_ref = 503    # varies depending on the reference
gwas_thresh = 5e-8  
single_snp_heidi_thresh = 0.01
multi_snps_heidi_thresh = 0.01
nsnps_thresh = 0   # the minimum number of instruments required for the GSMR analysis
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.03    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development 
gsmr_results_CMD = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 

#res
cat("The estimated effect of the exposure on outcome: ",exp(gsmr_results_CMD$bxy))
cat("95%L of bxy: ",exp(gsmr_results_CMD$bxy - 1.96 * gsmr_results_CMD$bxy_se))
cat("95%u of bxy: ",exp(gsmr_results_CMD$bxy + 1.96 * gsmr_results_CMD$bxy_se))
cat("P-value for bxy: ", gsmr_results_CMD$bxy_pval)
cat("Number of the SNPs used in the GSMR analysis: ", length(gsmr_results_CMD$used_index))
cat("Number of SNPs with missing estimates in the summary data: ", length(gsmr_results_CMD$na_snps))
cat("Number of non-significant SNPs: ", length(gsmr_results_CMD$weak_snps))
cat("Number of SNPs in high LD ( LD rsq >", ld_r2_thresh, "): ", length(gsmr_results_CMD$linkage_snps))
cat("Number of pleiotropic outliers: ", length(gsmr_results_CMD$pleio_snps))