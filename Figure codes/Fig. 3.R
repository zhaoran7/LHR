#---------------------------------------
# Fig. 3
#---------------------------------------
dat <- read_excel('/Source Data/Source Data Fig. 3.xlsx')

outcome_order <- unique(dat$Outcome)
plot_data$Outcome <- factor(dat$Outcome, levels = outcome_order)

ggplot(plot_data, aes(x = factor(height.cat), y = HR, fill = factor(lhr.cat))) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.75, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.2) +
  geom_text(aes(label = sig),
            position = position_dodge(0.8), vjust = -0.8, size = 2.8, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(-Outcome, ncol = 3) +
  scale_fill_manual(values = c("1" = "#F9BEB9", "2" = "#F2F2F2", "3" = "#AFC7EA"),
                    labels = c("Low LHR", "Medium LHR", "High LHR")) +
  labs(x = "Height Category", y = "Hazard Ratio (HR)", fill = "LHR") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "right"
  )

