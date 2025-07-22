#---------------------------------------
# Fig. 1c, 1d, and Supplementary Fig. 4b
#---------------------------------------
library(data.table);library(ggplot2);library(dplyr)

# height of our GWAS and GIANT 2022 as example
dat1 <- fread("/data/GWAS/height/height.gz") 
dat1 <- dat1 %>% filter(P < 5e-8) %>% select(SNP, EA, NEA, beta1 = BETA, p1 = P)

dat2 <- fread("/data/GWAS/height/HEIGHT_GIANT_2022.gz") 
dat2 <- dat2 %>% filter(P < 5e-8) %>% select(SNP, EA, NEA,beta2 = BETA, p2 = P)

dat <- inner_join(dat1, dat2, by = "SNP")
dat <- dat %>% filter((EA.x == EA.y & NEA.x == NEA.y) | (EA.x == NEA.y & NEA.x == EA.y)) %>%
  mutate(beta2_aligned = ifelse(EA.x == EA.y, beta2, -beta2))

ggplot(dat, aes(x = beta1, y = beta2_aligned)) +
  geom_point(alpha = 0.6, color = "#3167CD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
  xlim(-0.2, 0.25) + ylim(-0.2, 0.2) +
  labs(
    x = expression("SNP effect on (per SD) on Height in our GWAS"),
    y = expression("SNP effect on (per SD) on Height in GIANT (2022)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + annotate("text", x = -0.1, y = -0.15, label = "Genetic correlation (rg = 0.97)", size = 4, hjust = 0)
ggsave("height-leg.pdf",width = 5, height = 4.5)
