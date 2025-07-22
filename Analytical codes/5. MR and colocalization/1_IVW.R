#---------------------------------------
# Mendelian randomization (IVW method)
#---------------------------------------
pacman::p_load(data.table, dplyr, stringi, TwoSampleMR)

# gwas files
exposure_list <- list(lhr = "LHR.gz",height="height.gz")
outcome_dir <- "/data/GWAS/CMD/"
outcome_files <- list.files(outcome_dir, pattern = "*.gz", full.names = TRUE)
all_results <- list()

# loop
for (x in names(exposure_list)) {
  print(paste('processing',x))
  dat.X.raw <- fread(paste0("/data/GWAS/height/", exposure_list[[x]]))
  # excluding
  dat.X.raw$P <- as.numeric(dat.X.raw$P); dat.X.sig <- dat.X.raw %>% filter(P <= 5e-8)
  dat.X.sig <- dat.X.sig %>% filter(!(CHR == 6 & POS >= 28477797 & POS <= 33448354))
  dat.X.sig <- dat.X.sig %>% mutate(F_stat = (BETA^2) / (SE^2)) %>% filter(F_stat >= 10)
  # independent SNPs after conditional analysis
  dat.X.top <- fread(paste0("/data/GWAS/height/", x, ".jma.cojo"))
  dat.X <- data.frame(dat.X.sig[SNP %in% dat.X.top$SNP, ]) %>%
    format_data(type = 'exposure',snp_col = 'SNP', effect_allele_col = 'EA', other_allele_col = 'NEA',
                beta_col = 'BETA', se_col = 'SE', pval_col = 'P',samplesize_col = 'N', eaf_col = 'EAF')
  
  for (yfile in outcome_files) {
    y <- tools::file_path_sans_ext(basename(yfile))
    print(paste('processing',y))
    dat.Y.raw <- fread(yfile)
    dat.Y <- data.frame(dat.Y.raw %>% merge(dat.X, by = "SNP")) %>%
      format_data(type = 'outcome',snp_col = 'SNP', effect_allele_col = 'EA', other_allele_col = 'NEA',
                  beta_col = 'BETA', se_col = 'SE', pval_col = 'P',samplesize_col = 'N', eaf_col = 'EAF')
    
    dat.CVD <- harmonise_data(dat.X, dat.Y, action = 1)
    res <- mr(dat.CVD) %>% generate_odds_ratios()
    ivw_res <- res %>% filter(method == "Inverse variance weighted") %>% mutate(exposure = x, outcome = y)
    all_results[[paste(x, y, sep = "_")]] <- ivw_res
  }
}
final_df <- bind_rows(all_results)
fwrite(final_df, "/res/mr_ivw.csv")