pacman::p_load(clusterProfiler, org.Hs.eg.db, openxlsx, dplyr, stringr)

gene_dir <- "/data/magma/shared_genes"
out_dir <- "/data/magma"

#---------------------------------------
# step 1: GO and KEGG for Height and LHR
#---------------------------------------
traits <- c("Height", "LHR")
trait_files <- c(
  Height = file.path(out_dir, "Source Data Supplementary Fig. 14a.xlsx"),
  LHR    = file.path(out_dir, "Source Data Supplementary Fig. 14b.xlsx")
)

for (trait in traits) {
  genes <- scan(file.path(gene_dir, paste0(trait, ".sig.genes.txt")), what = "character")
  wb <- createWorkbook()
  
  for (ont in c("BP", "CC", "MF")) {
    ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = ont, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    pAdjustMethod = "BH", readable = TRUE)
    addWorksheet(wb, paste0("GO_", ont))
    writeData(wb, paste0("GO_", ont), as.data.frame(ego))
  }
  
  ekegg <- enrichKEGG(genes, organism = "hsa", pvalueCutoff = 0.05)
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  addWorksheet(wb, "KEGG")
  writeData(wb, "KEGG", as.data.frame(ekegg))
  
  saveWorkbook(wb, trait_files[[trait]], overwrite = TRUE)
}

#---------------------------------------
# step 1: GO and KEGG for shared genes
#---------------------------------------
setwd(gene_dir)
file_list <- list.files(pattern = "\\.shared\\.genes\\.txt$")
go <- list(); kegg <- list()

for (file in file_list) {
  genes <- scan(file, what = "character")
  if (length(genes) < 10) next
  ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH", readable = TRUE)
  go_df <- as.data.frame(ego); if (nrow(go_df) > 0) go[[file]] <- mutate(go_df, File = file)
  
  ekegg <- enrichKEGG(genes, organism = "hsa", pvalueCutoff = 0.05)
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(ekegg); if (nrow(kegg_df) > 0) kegg[[file]] <- mutate(kegg_df, File = file)
}

go_all <- bind_rows(go)
kegg_all <- bind_rows(kegg)

wb2 <- createWorkbook()
addWorksheet(wb2, "GO_BP")
addWorksheet(wb2, "KEGG")
writeData(wb2, "GO_BP", go_all)
writeData(wb2, "KEGG", kegg_all)
saveWorkbook(wb2, file.path(out_dir, "Source Data Fig. 7.xlsx"), overwrite = TRUE)
