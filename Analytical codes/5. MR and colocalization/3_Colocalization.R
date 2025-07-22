#---------------------------------------
# colocalization
#---------------------------------------
pacman::p_load(coloc,data.table,ieugwasr,dplyr)
# Define window
GWAS_WINDOW <- 500000 

# LHR and T2D colocalization
# Load X
dat.X.raw <- fread("/data/GWAS/height/LHR.gz") 
# Load Y
dat.Y.raw <- fread("/data/GWAS/CMD/T2D.gz",fill=T)  

# ====== X significant SNPs
dat.X.sig <- dat.X.raw[which(dat.X.raw$P < 5e-8), ]
## Use ieugwasr for LD clumping
dat.X.sigsnp <- ieugwasr::ld_clump_local(
  data.frame(rsid = dat.X.sig$SNP, pval = dat.X.sig$P),clump_kb = 1000, clump_r2 = 0.001, clump_p = 1,
  bfile = "/data/1kG/EUR",plink_bin = "/software/plink_mac_20241022/plink")
dat.X.clump <- dat.X.sig[SNP %in% dat.X.sigsnp$rsid, ]
dat.X.clump[, `:=`(start = POS - GWAS_WINDOW, end = POS + GWAS_WINDOW)]

##====== Y significant SNPs
dat.Y.sig <- dat.Y.raw[which(dat.Y.raw$P < 5e-8), ]
## Use ieugwasr for LD clumping
dat.Y.sigsnp <- ieugwasr::ld_clump_local(
  dat = data.frame(rsid = dat.Y.sig$SNP, pval = dat.Y.sig$P), clump_kb = 1000, clump_r2 = 0.1, clump_p = 1,
  bfile = "/data/1kG/EUR",plink_bin = "/software/plink_mac_20241022/plink")
dat.Y.clump <- dat.Y.sig[SNP %in% dat.Y.sigsnp$rsid, ]
dat.Y.clump[, `:=`(start = POS - GWAS_WINDOW, end = POS + GWAS_WINDOW)]

#SNPs for coloc
glist <- rbind(dat.X.clump[, .(CHR, start, end, SNP)],dat.Y.clump[, .(CHR, start, end, SNP)])
#result list
coloc_results <- list()

for (i in 1:nrow(glist)) {
  curr_snp <- glist$SNP[i];curr_chr <- glist$CHR[i];curr_start <- glist$start[i];curr_end <- glist$end[i]
  tryCatch({
    ## Extract region for height and AF data
    X_curr_region <- dat.X.raw[CHR == curr_chr & POS >= curr_start & POS <= curr_end, ]
    Y_curr_region <- dat.Y.raw[CHR == curr_chr & POS >= curr_start & POS <= curr_end, ]
    ## Skip if no SNPs in either dataset
    if (nrow(X_curr_region) == 0 | nrow(Y_curr_region) == 0) {
      message(paste("Skipping SNP", curr_snp, "due to no SNPs in one or both datasets"))
      next
    }
    ## Merge and deduplicate
    data_merged <- merge(X_curr_region,Y_curr_region,by = c("CHR", "POS"),suffixes = c("_X", "_Y"))
    data_merged <- data_merged[!duplicated(data_merged$SNP_X), ]
    ## Align effect alleles
    data_merged <- data_merged %>% dplyr::filter((EA_X == EA_Y & NEA_X == NEA_Y) | (EA_X == NEA_Y & NEA_X == EA_Y))
    data_merged <- data_merged %>% dplyr::mutate(BETA_Y = ifelse(EA_X == EA_Y, BETA_Y, -BETA_Y))
    ## Calculate variance
    data_merged$VAR_X <- data_merged$SE_X^2
    data_merged$VAR_Y <- data_merged$SE_Y^2
    data_merged <- data_merged[data_merged$VAR_X != 0 & data_merged$VAR_Y != 0, ]
    
    ## Skip if no SNPs after filtering
    if (nrow(data_merged) == 0) {
      message(paste("Skipping SNP", curr_snp, "due to no SNPs after filtering"))
      next
    } 
    ## Split and format data for coloc
    X_coloc <- data_merged[, .(beta = BETA_X, varbeta = VAR_X, snp = SNP_X, MAF = EAF_X, N = N_X)]
    Y_coloc <- data_merged[, .(beta = BETA_Y, varbeta = VAR_Y, snp = SNP_Y, MAF = EAF_Y, N = N_Y)]
    X_coloc <- as.list(X_coloc)
    Y_coloc <- as.list(Y_coloc)
    ## Declare phenotype types
    X_coloc$type <- "quant"  
    Y_coloc$type <- "cc"    
    ## Run coloc analysis
    res <- coloc.abf(X_coloc, Y_coloc, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    coloc_results[[curr_snp]] <- res$summary
  }, error = function(e) {
    message(paste("Error in SNP", curr_snp, ":", e$message))
  })
}
coloc_results_df <- do.call(rbind, lapply(names(coloc_results), function(snp) {
  data.frame(SNP = snp, t(coloc_results[[snp]]))
}))
fwrite(coloc_results_df,'coloc.csv',row.names=F)
