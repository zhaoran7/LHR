#---------------------------------------
# Fig. 2 and Supplementary Fig. 9,10,12
#---------------------------------------
library(readxl)
library(dplyr)
library(forestploter)
forest <- read_excel('/Source Data/Source Data Fig. 2.xlsx',sheet=1)

forest$CMD <- ifelse(is.na(forest$CMD),"",forest$CMD);forest <- forest %>% add_column(X = NA, .after = "CMD")
forest$X <- paste(rep("   ", 15), collapse = "  ");colnames(forest)[3] <- " "

tm <- forest_theme(base_size = 8,
                   ci_pch = 16,
                   ci_col = "#00539f",
                   ci_fill = "#00539f",
                   refline_lty = "solid",
                   refline_lwd = 0.5,
                   legend_value = "exposure",
                   arrow_type = "closed",
                   xaxis_lwd = 0.6,
                   xaxis_cex = 0.8,
                   title_just = "center",
                   core=list(bg_params=list(fill = c("#DDF1FF","#DDF1FF","white","white"))))
forest.plot <- forest((forest[,c(2,3,5,9)]),
                      est=list(forest$HR),
                      lower=list(forest$Lower),
                      upper=list(forest$Upper),
                      sizes = 0.6,
                      ref_line = 1,
                      xlim = c(0.7, 1.5),
                      ticks_at = c(0.8,1,1.2,1.4),
                      ci_column=c(2),
                      nudge_y = 0,
                      theme=tm)
print(forest.plot)
