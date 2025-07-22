#---------------------------------------
# Fig. 5a
#---------------------------------------
pacman::p_load(readxl,dplyr,tidyr,circlize,ComplexHeatmap)

dat <- read_excel('/Source/Source Data Fig. 5.xlsx', sheet = 1)
dat <- dat %>% mutate(Trait = paste0(exposure, "-", outcome))
dat <- dat %>%
  mutate(snp = factor(snp, levels = unique(snp)))
dat_mat <- dat %>%
  dplyr::select(Trait, snp, PP.H4) %>%
  pivot_wider(names_from = Trait, values_from = PP.H4) %>%
  column_to_rownames("snp") %>%
  as.matrix()
dat_mat <- dat_mat[rowSums(is.na(dat_mat)) < ncol(dat_mat), ]
dat_mat <- dat_mat[, colSums(is.na(dat_mat)) < nrow(dat_mat)]
dat_mat <- dat_mat[rev(rownames(dat_mat)),]

green_pink <- colorRamp2(breaks = c(0.9, 1.0), colors = c("#ebbfc2", "#e28187"))

# circos param
circos.clear()
circos.par(start.degree = 90, gap.degree = 30, clock.wise = FALSE,
           canvas.xlim = c(-1, 1), canvas.ylim = c(-0.3, 0.3),
           track.margin = c(0, 0.01),
           cell.padding = c(0, 0, 0, 0))

# heatmap
circos.heatmap(dat_mat,
               cluster = FALSE,
               bg.border = "black", bg.lwd = 0.5, na.col = "white",
               cell.border = "grey", cell.lwd = 0.5,
               rownames.side = "outside", rownames.cex = 0.7,
               col = green_pink,
               track.height = 0.2)

#  Trait 
circos.track(track.index = get.current.track.index(),
             bg.border = NA,
             panel.fun = function(x, y) {
               cn <- rev(colnames(dat_mat))
               n <- length(cn)
               cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
               y_coords <- seq(CELL_META$cell.ylim[1] + cell_height / 2,
                               CELL_META$cell.ylim[2] - cell_height / 2,
                               length.out = n)
               for (i in 1:n) {
                 circos.lines(
                   c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
                   c(y_coords[i], y_coords[i]), col = "black", lwd = 1)
               }
               circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
                           y_coords, labels = cn, cex = 0.6, adj = c(0, 0.5), facing = "inside")
             })

# legend
lgd <- Legend(col_fun = green_pink, title = "PP.H4", at = c(0.9, 0.95, 1))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")

#---------------------------------------
# Fig. 5b-f and Supplementary Fig. 2
#---------------------------------------
pacman::p_load(loucszoomr,LDlinkR,EnsDb.Hsapiens.v86,data.table)

# rs66922415 as example
data.X.raw <- fread("/data/GWAS/height/height.gz")
data.M.raw <- fread("/data/GWAS/height/LHR.gz")
data.Y.raw <- fread("/data/GWAS/CVD/T2D.gz")

WINDOW <- 500000
top_snp <-"rs66922415"

top_snp_info <- data.X.raw[SNP == top_snp, .(CHR, POS)]
top_snp_chr <- top_snp_info$CHR
top_snp_pos <- top_snp_info$POS

region_X <- data.X.raw[CHR == top_snp_chr & POS >= (top_snp_pos - WINDOW) & POS <= (top_snp_pos + WINDOW), ]
region_M <- data.M.raw[CHR == top_snp_chr & POS >= (top_snp_pos - WINDOW) & POS <= (top_snp_pos + WINDOW), ]
region_Y <- data.Y.raw[CHR == top_snp_chr & POS >= (top_snp_pos - WINDOW) & POS <= (top_snp_pos + WINDOW), ]

common_snps <- intersect(region_X$SNP, intersect(region_M$SNP, region_Y$SNP))
writeLines(common_snps, "common_snps.txt")

#LD R2
#plink --bfile /Volumes/PS2000/Height2CVD/data/magma/ref/g1000_eur/g1000_eur \
#--ld-snp rs66922415 \
#--extract common_snps.txt \
#--ld-window 99999 \
#--ld-window-kb 1000 \
#--ld-window-r2 0 \
#--out r2_values \
#--r2

region_X <- region_X[SNP %in% common_snps, ]
region_M <- region_M[SNP %in% common_snps, ]
region_Y <- region_Y[SNP %in% common_snps, ]
r2_values <- read.table("r2_values.ld", header = TRUE)
r2_values <- r2_values[, (ncol(r2_values)-1):ncol(r2_values)]
colnames(r2_values) <- c("rsid", "R2")

rename_and_merge <- function(data) {
  setnames(data, old = c("SNP", "CHR", "POS", "EA", "NEA", "EAF", "N", "BETA", "SE", "P"), 
           new = c("rsid", "chr", "pos", "ea", "nea", "eaf", "n", "beta", "se", "p"), 
           skip_absent = TRUE)
  merge(data, r2_values, by.x = "rsid", by.y = "rsid", all = FALSE)
}

region_X <- rename_and_merge(region_X)
region_M <- rename_and_merge(region_M)
region_Y <- rename_and_merge(region_Y)

# draw region_X, region_M, region_Y separately
loc <- locus(data = region_X, index_snp = top_snp, flank = 25e4, ens_db = "EnsDb.Hsapiens.v86", LD = "R2")
color_mapping <- c("grey95","#87ceeb","#006401","#ffa501","#ff0000")
loc$data$bg <- ifelse(is.na(loc$data$ld),"#e9e9e9",color_mapping[findInterval(loc$data$ld, seq(0, 1, 0.2), rightmost.closed = TRUE)])
loc$data$col <- "grey95"
index_snp <- which(loc$data$rsid == loc$index_snp)
loc$data$bg[index_snp] <- "#a01ff0"
loc$data$cex <- 1

region_X <- locus_plot(loc,labels = top_snp,border = TRUE, cex.axis =0.8,cex.lab = 1,heights = c(1, 0.5),
           cex.text = 0.4,italics = TRUE,text_pos = "top",exon_border = NA,xticks = "bottom",filter_gene_name = c('MIA3'))
ggsave('region_X.pdf',region_X, width =4,height=3)