pacman::p_load(dplyr, readxl, ggplot2, patchwork, circlize, ComplexHeatmap, RColorBrewer)
file_path <- "/Volumes/PS2000/Height2CVD/manuscript/NC/Source Data Fig. 6.xlsx"

#---------------------------------------
# Fig. 6a
#---------------------------------------
df <- read_excel(file_path, sheet = 1)

df_height <- df %>%
  filter(Trait == "height") %>%
  arrange(desc(input_list_combined_log10p)) %>%
  distinct(General_cell_type, .keep_all = TRUE) %>%
  slice_head(n = 15) %>%
  mutate(General_cell_type = factor(General_cell_type, levels = General_cell_type))

df_lhr <- df %>%
  filter(Trait == "LHR") %>%
  arrange(desc(input_list_combined_log10p)) %>%
  distinct(General_cell_type, .keep_all = TRUE) %>%
  slice_head(n = 15) %>%
  mutate(General_cell_type = factor(General_cell_type, levels = General_cell_type))

p1 <- ggplot(df_height, aes(x = General_cell_type, y = input_list_combined_log10p)) +
  geom_bar(stat = "identity", fill = "#3B9DFF") +
  labs(title = "Height", x = "Cell Type", y = "-log10(adjP)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

p2 <- ggplot(df_lhr, aes(x = General_cell_type, y = input_list_combined_log10p)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  labs(title = "LHR", x = "Cell Type", y = "-log10(adjP)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

(p1 / p2) + plot_layout(heights = c(1, 1)) %>%
  ggsave("Fig. 6a.pdf", ., width = 7, height = 9)

#---------------------------------------
# Fig. 6b
#---------------------------------------
df <- read_excel(file_path, sheet = 2)

df <- df %>%
  filter(p < 0.05) %>%
  distinct(Trait, GeneSet, N_genes, p) %>%
  group_by(Trait, GeneSet) %>%
  summarise(p_bonf = min(p) * n(), .groups = "drop") %>%
  mutate(p_bonf = ifelse(p_bonf > 1, 1, p_bonf),
         log10P = -log10(p_bonf)) %>%
  select(Trait, GeneSet, log10P) %>%
  group_by(Trait) %>%
  slice_max(order_by = log10P, n = 10, with_ties = FALSE) %>%
  ungroup()

trait_colors <- setNames(
  rep(brewer.pal(12, "Set3"), length.out = length(unique(df$Trait))),
  unique(df$Trait)
)
tissue_colors <- setNames(rep("grey70", length(unique(df$GeneSet))), unique(df$GeneSet))
grid_col <- c(tissue_colors, trait_colors)
link_col <- trait_colors[df$Trait]

pdf("Fig. 6b.pdf", width = 7, height = 7)
plot.new()
circos.clear()
chordDiagramFromDataFrame(df, grid.col = grid_col, col = link_col,
                          annotationTrack = "grid",
                          annotationTrackHeight = c(0.05, 0.03),
                          transparency = 0.2, directional = 0)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(mean(get.cell.meta.data("xlim")),
              get.cell.meta.data("ylim")[1] + 1.2,
              get.cell.meta.data("sector.index"),
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

legend_weights <- c(1.3, 10, 20)
leg <- Legend(title = "-log10P",
              labels = as.character(legend_weights),
              legend_gp = gpar(lwd = legend_weights / 2 , col = "black"),
              type = "lines")
draw(packLegend(leg), x = unit(1, "npc") - unit(10, "mm"), just = "right")
dev.off()

#---------------------------------------
# Fig. 6c
#---------------------------------------
df <- read_excel(file_path, sheet = 3)

top15_height <- df %>%
  filter(Trait == "height") %>%
  arrange(desc(input_list_combined_log10p)) %>%
  slice_head(n = 15) %>%
  mutate(General_cell_type = factor(General_cell_type, levels = rev(unique(General_cell_type))))

top15_lhr <- df %>%
  filter(Trait == "LHR") %>%
  arrange(desc(input_list_combined_log10p)) %>%
  slice_head(n = 15) %>%
  mutate(General_cell_type = factor(General_cell_type, levels = rev(unique(General_cell_type))))

p1 <- ggplot(top15_height, aes(x = General_cell_type, y = input_list_combined_log10p)) +
  geom_bar(stat = "identity", fill = "#3B9DFF") +
  labs(title = "Height", x = "Tissue", y = "-log10(adjP)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

p2 <- ggplot(top15_lhr, aes(x = General_cell_type, y = input_list_combined_log10p)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  labs(title = "LHR", x = "Tissue", y = "-log10(adjP)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

(p1 / p2) + plot_layout(heights = c(1, 1)) %>%
  ggsave("Fig. 6c.pdf", ., width = 7, height = 9)

#---------------------------------------
# Fig. 6d
#---------------------------------------
df <- read_excel(file_path, sheet = 4)

df_chord <- df %>%
  filter(input_list_combined_log10p >= 1.3) %>%
  group_by(Trait) %>%
  slice_max(order_by = input_list_combined_log10p, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  rename(Cell = General_cell_type) %>%
  mutate(weight = input_list_combined_log10p) %>%
  select(Trait, Cell, weight)

trait_colors <- setNames(
  rep(brewer.pal(12, "Set3")[-2], length.out = length(unique(df_chord$Trait))),
  unique(df_chord$Trait)
)
cell_colors <- setNames(rep("grey70", length(unique(df_chord$Cell))), unique(df_chord$Cell))
grid_col <- c(trait_colors, cell_colors)
link_col <- trait_colors[df_chord$Trait]

pdf("Fig. 6d.pdf", width = 7, height = 7)
plot.new()
circos.clear()
chordDiagramFromDataFrame(df_chord, grid.col = grid_col, col = link_col,
                          annotationTrack = "grid",
                          annotationTrackHeight = c(0.05, 0.03),
                          transparency = 0.4, directional = 0)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(mean(get.cell.meta.data("xlim")),
              get.cell.meta.data("ylim")[1] + 1.2,
              get.cell.meta.data("sector.index"),
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

legend_weights <- c(1.3, 2, 5)
leg <- Legend(title = "-log10(combined P)",
              labels = as.character(legend_weights),
              legend_gp = gpar(lwd = legend_weights, col = "black"),
              type = "lines")
draw(packLegend(leg), x = unit(1, "npc") - unit(10, "mm"), just = "right")
dev.off()
