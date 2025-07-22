#---------------------------------------
# Supplementary Fig. 14
#---------------------------------------
pacman::p_load(dplyr, ggplot2, openxlsx, patchwork)

plot_list <- list()

for (trait in c("Height", "LHR")) {
  file <- sprintf("/Source Data/Source Data Supplementary Fig. 14%s.xlsx",ifelse(trait == "Height", "a", "b"))
  sheets <- getSheetNames(file)
  
  dfs <- lapply(sheets, function(s) read.xlsx(file, sheet = s))
  names(dfs) <- sheets
  for (i in names(dfs)) dfs[[i]]$Category <- sub("GO_", "", i)
  
  all_top <- bind_rows(dfs) %>%
    arrange(factor(Category, levels = c("GO BP", "GO CC", "GO MF", "KEGG")), desc(Count)) %>%
    group_by(Category) %>% slice_head(n = 10) %>% ungroup() %>%
    mutate(log10p = -log10(p.adjust), Description = factor(Description, levels = rev(unique(Description)))) %>%
    mutate(yid = as.numeric(Description))
  
  cat_colors <- c("GO BP"="#ea739c", "GO CC"="#eeba54", "GO MF"="#b4d88c", "KEGG"="#a0c8e8")
  cat_ranges <- all_top %>% group_by(Category) %>%
    summarise(y_min = min(yid), y_max = max(yid), .groups = "drop")
  plot_titles <- c("Height" = "a", "LHR" = "b")
  p <- ggplot(all_top, aes(y = Description)) +
    geom_rect(data = cat_ranges, aes(xmin=-2.5, xmax=-1.8, ymin=y_min-0.4, ymax=y_max+0.4, fill=Category),
              color=NA, inherit.aes=FALSE, show.legend=FALSE) +
    geom_point(aes(x=-1, size=Count, fill=Category), shape=21, color="white", stroke=1.2) +
    geom_text(aes(x=-1, label=Count), size=2.5, color="black") +
    geom_col(aes(x=log10p, fill=Category), width=0.65, show.legend=FALSE) +
    geom_text(aes(x=0.22, label=Description), hjust=0, size=2.5, color="black") +
    scale_fill_manual(values = cat_colors) +
    scale_size_continuous(range = c(4, 6)) +
    scale_x_continuous(expand = expansion(mult = c(0.22, 0.02)),
                       name = expression(-log[10](adjusted-italic(P)))) +
    ggtitle(plot_titles[[trait]]) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.6),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(10, 15, 10, 30))
  
  plot_list[[trait]] <- p
}
#save
combined_plot <- plot_list$Height | plot_list$LHR + plot_layout(nrow = 1)
ggsave("Supplementary Fig. 14.pdf", plot = combined_plot, width = 10, height = 16)
