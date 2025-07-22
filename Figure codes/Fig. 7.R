#---------------------------------------
# Fig. 7
#---------------------------------------
pacman::p_load(dplyr, ggplot2, stringr, openxlsx)

# data
data_file <- "/Source Data/Source Data Fig. 7.xlsx"
go_all <- read.xlsx(data_file, sheet = "Source Data Fig. 7a")
kegg_all <- read.xlsx(data_file, sheet = "Source Data Fig. 7b")

# GO
go_plot <- go_all %>%
  mutate(mLog10P = -log10(p.adjust)) %>%
  filter(Trait %in% c("Height-AF", "Height-CAD", "Height-HT", "Height-T2D",
                         "LHR-CAD", "LHR-HT", "LHR-MI", "LHR-T2D")) %>%
  filter(!str_detect(Description, "MHC")) %>%
  group_by(Trait) %>%
  slice_max(mLog10P, n = 10) %>%
  ungroup()

p_go <- ggplot(go_plot, aes(x = Trait, y = Description)) +
  geom_point(aes(size = Count, color = mLog10P)) +
  scale_color_gradient(low = "#ebbfc2", high = "#904550") +
  scale_size_continuous(range = c(2, 6)) +
  labs(y = NULL, x = NULL, color = "-log10(P adj)", size = "Count") +
  ggtitle("a") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0, size = 14),
        legend.position = "right")

# KEGG
kegg_plot <- kegg_all %>%
  mutate(mLog10P = -log10(p.adjust)) %>%
  filter(Trait %in% c("Height-AF", "Height-CAD", "Height-HT", "Height-T2D",
                         "LHR-CAD", "LHR-HT", "LHR-MI", "LHR-T2D")) %>%
  filter(category != "Human Diseases") %>%
  group_by(Trait) %>%
  slice_max(mLog10P, n = 15) %>%
  ungroup()

p_kegg <- ggplot(kegg_plot, aes(x = Trait, y = Description)) +
  geom_point(aes(size = Count, color = mLog10P)) +
  scale_color_gradient(low = "#8FB4BE", high = "#0B425E") +
  scale_size_continuous(range = c(2, 6)) +
  labs(y = NULL, x = NULL, color = "-log10(P adj)", size = "Count") +
  ggtitle("b") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0, size = 14),
        legend.position = "right")

combined_plot <- p_go | p_kegg
ggsave("Fig. 7.pdf", combined_plot, width = 16, height = 8)
