#---------------------------------------
# Supplementary Fig. 8
#---------------------------------------
pacman::p_load(corrplot, ggplot2, dplyr)

dat <- readRDS("/Volumes/PS2000/Height2CVD/data/height.RDS")

# Supplementary Fig. 8a
cor_mat <- cor(dat[, c("height", "leg", "LHR")], use = "complete.obs")
corrplot(cor_mat, method = "circle", type = "lower", addCoef.col = "black", tl.col = "black")

# Supplementary Fig. 8b
dat %>%
  mutate(z_height = scale(height)[, 1],
         z_lhr = scale(lhr.res)[, 1]) %>%
  ggplot(aes(x = z_height, y = z_lhr)) +
  geom_point(alpha = 0.3, size = 0.6, color = "#4682B4") +
  labs(x = "Standardized Height (Z-score)",
       y = "Standardized LHR Residuals (Z-score)") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-10, 10)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.22)
  )