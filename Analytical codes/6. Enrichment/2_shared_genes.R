#---------------------------------------
# shared genes
#---------------------------------------
pacman::p_load(dplyr, AnnotationDbi, org.Hs.eg.db)

shared_dir <- "/data/magma/shared_genes"
out_dir <- "/data/magma/output"
dir.create(shared_dir, showWarnings = FALSE)

#---------------------------------------
# step 1: convert ENTREZ → SYMBOL
#---------------------------------------
gene_files <- list.files(shared_dir, pattern = "\\.shared\\.genes\\.txt$", full.names = TRUE)

for (file in gene_files) {
  cat("Processing:", file, "\n")
  if (file.info(file)$size == 0) {
    cat("⚠️  Skipped (empty):", file, "\n\n")
    next
  }
  ids <- scan(file, what = "", quiet = TRUE)
  if (length(ids) == 0) {
    cat("⚠️  Skipped (no IDs):", file, "\n\n")
    next
  }
  
  sym <- AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns = "SYMBOL", keytype = "ENTREZID") %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(SYMBOL) %>%
    pull(SYMBOL)
  
  out <- sub("\\.shared\\.genes\\.txt$", ".symbol.txt", file)
  writeLines(sym, out)
  cat("Saved:", out, "(n =", length(sym), "genes)\n\n")
}

cat("✅ All SYMBOL conversion done! 🚀\n")

#---------------------------------------
# step 2: find shared genes
#---------------------------------------
bonf_cutoff <- 0.05
trait_list <- c("AF","VTE","CM","AA","T2D","CAD","HF","CKD","HS","IS","PAD","AS","TIA","CAD","HT","MI")

get_sig_genes <- function(file, cutoff) {
  dat <- read.table(file, header = TRUE)
  dat$BONF <- p.adjust(dat$P, method = "bonferroni")
  dat$GENE[dat$BONF < cutoff]
}

# height and LHR
height_genes <- get_sig_genes(file.path(out_dir, "Height.genes.out"), bonf_cutoff)
lhr_genes <- get_sig_genes(file.path(out_dir, "LHR.genes.out"), bonf_cutoff)

writeLines(height_genes, file.path(shared_dir, "Height.sig.genes.txt"))
writeLines(lhr_genes, file.path(shared_dir, "LHR.sig.genes.txt"))

# CMD traits
for (trait in trait_list) {
  cat("Processing:", trait, "\n")
  trait_file <- file.path(out_dir, paste0(trait, ".genes.out"))
  if (!file.exists(trait_file)) {
    cat("⚠️  File not found:", trait_file, "\n\n")
    next
  }
  
  trait_genes <- get_sig_genes(trait_file, bonf_cutoff)
  writeLines(trait_genes, file.path(shared_dir, paste0(trait, ".sig.genes.txt")))
  
  writeLines(intersect(height_genes, trait_genes),
             file.path(shared_dir, paste0("Height_vs_", trait, ".shared.genes.txt")))
  writeLines(intersect(lhr_genes, trait_genes),
             file.path(shared_dir, paste0("LHR_vs_", trait, ".shared.genes.txt")))
}
