#---------------------------------------
# Fig. 1b
#---------------------------------------

library(corrplot);library(readxl)
dat <- read_excel('/Source Data/Source Data Fig. 1b.xlsx')
 dat[c("Trait 1", "Trait 2")] <- lapply(dat[c("Trait 1", "Trait 2")], function(x) gsub("height", "Height", gsub("lhr.res", "LHR", x)))
 traits <- union(dat$'Trait 1', dat$'Trait 2')
mat <- matrix(NA, nrow = length(traits), ncol = length(traits),dimnames = list(traits, traits))
for (i in seq_len(nrow(dat))) {mat[dat$'Trait 2'[i], dat$'Trait 1'[i]] <- dat$rg[i]}
mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)];diag(mat) <- 1
order <- c('Height','LHR','AF','VTE','CM','AA','HF','CKD','HS','IS','PAD','AS','T2D','TIA','CAD','HT','MI')
mat <- mat[order, order]

dat$p_fdr <- p.adjust(dat$p, method = "fdr")
p_mat <- matrix(NA, nrow = length(traits), ncol = length(traits), dimnames = list(traits, traits))
for (i in seq_len(nrow(dat))) {p_mat[dat$'Trait 2'[i], dat$'Trait 1'[i]] <- dat$p_fdr[i]}
p_mat[upper.tri(p_mat)] <- t(p_mat)[upper.tri(p_mat)];diag(p_mat) <- 0
p_mat <- p_mat[order, order]

pdf("Fig. 1b.pdf", width = 8, height = 8)
corrplot(mat, method = "square", order = "original",type = "lower", number.cex = 0.5,
         p.mat = p_mat, sig.level = c(0.001,0.01,0.05), insig = 'label_sig',
         tl.col = "black", tl.cex = 0.7, pch.cex =0.5,tl.srt = 45, tl.pos = "lt")
corrplot(mat, method = "number",order = "original", type = "upper", 
         tl.col = "n", tl.cex = 0.6,number.cex = 0.5, tl.pos = "n",add = TRUE )
dev.off()