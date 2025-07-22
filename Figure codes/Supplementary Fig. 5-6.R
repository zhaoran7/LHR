#---------------------------------------
# Supplementary Fig. 5-6
#---------------------------------------
pacman::p_load(dplyr, stringr, readr, tidyr, ggplot2, readxl)

# Read annotation mapping
# download: https://console.cloud.google.com/storage/browser/_details/broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores/biorxiv/GTEx_1000Gv3.tgz
annot_map <- read.delim("/Users/zhaor/Downloads/LDSCORE_LDSC_SEG_ldscores_biorxiv_GTEx.ldcts", header = FALSE)
colnames(annot_map) <- c("TrueName", "AnnotPath")
annot_map <- annot_map %>%
  mutate(
    Annot = str_extract(AnnotPath, "GTEx\\.\\d+"),
    TrueName = str_replace_all(TrueName, "_", " ") %>%
      str_replace_all("\\(.*?\\)", "") %>%
      str_trim()
  )

# Read partitioned rg results
partitioned_ldsc <- read_xlsx('/Source Data/Source Data Supplementary Fig. 5.xlsx') %>%
  mutate(
    FDR_p = p.adjust(p, method = "fdr"),
    Trait = paste0(Trait2, "-", Trait1),
    LOG10P = -log10(FDR_p)
  ) %>%
  left_join(annot_map, by = "Annot") %>%
  mutate(Annot = coalesce(TrueName, Annot)) %>%
  select(-TrueName)

# Set factor orders
disease_order <- c('AF','VTE','CM','AA','HF','CKD','HS','IS','PAD','AS','T2D','TIA','CAD','HT','MI')
trait_order <- c(paste0("Height-", disease_order), paste0("LHR-", disease_order))
annot_order <- partitioned_ldsc %>% distinct(Annot) %>% arrange(tolower(Annot)) %>% pull(Annot)
annot_order <- c(setdiff(annot_order, "Total"), "Total")

partitioned_ldsc <- partitioned_ldsc %>%
  filter(p < 0.05) %>%
  mutate(
    Trait = factor(Trait, levels = trait_order),
    Annot = factor(Annot, levels = annot_order)
  )

# Heatmap of rg
ggplot(partitioned_ldsc %>%
         filter(p < 0.05) %>%
         complete(Trait = factor(Trait, levels = trait_order),
                  Annot = factor(Annot, levels = rev(annot_order))),
       aes(x = Trait, y = Annot, fill = rg)) +
  geom_tile(color = "grey80", linewidth = 0.2) +
  scale_fill_gradientn(
    name = "rg", limits = c(-0.33, 0.33),
    breaks = seq(-0.3, 0.3, 0.1),
    colours = c("#0B425E", "#8FB4BE", "white", "#E28187", "#904550"),
    values = scales::rescale(c(-0.33, -0.2, 0, 0.2, 0.33)),
    na.value = "white"
  ) +
  geom_text(aes(label = ifelse(is.na(rg), "", sprintf("%.2f", rg)),
                color = ifelse(abs(rg) >= 0.2, "white", "black")), size = 2.5) +
  scale_color_identity() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5, fill = NA),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.key = element_blank()
  ) +
  labs(x = NULL, y = NULL)

# Barplot of h2
h2_bar <- partitioned_ldsc %>%
  filter(grepl("^LHR-|^Height-", Trait)) %>%
  mutate(Trait2 = ifelse(grepl("^Height", Trait), "Height", "LHR")) %>%
  mutate(Trait2 = factor(Trait2, levels = c("LHR", "Height"))) %>%
  group_by(Annot, Trait2) %>%
  summarise(h2_1 = mean(h2_1, na.rm = TRUE), .groups = "drop") %>%
  mutate(Annot = factor(Annot, levels = rev(annot_order)))

ggplot(h2_bar, aes(x = h2_1, y = Annot, fill = Trait2)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Height" = "#1DBDC6", "LHR" = "#F6C490")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(x = "h2 (heritability)", fill = NULL)

