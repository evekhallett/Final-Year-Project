library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggpubr)
df <- read.csv('/Users/evehallett/Documents/dissertation/to_submit/Figure5/files/goodman_nt3.csv')
# Define the order of bases
base_order <- c('A', 'T', 'C', 'G')
plots <- list()
for (position in 2:11) {
  df[[paste0("Position_", position)]] <- factor(df[[paste0("Position_", position)]], levels = base_order)
  base_plot <- ggplot(df, aes_string(x = paste0("Position_", position), y = "Prot.FCC", fill = paste0("Position_", position))) +
    geom_violin() +
    labs(x = paste("Codon Position", position), y = "Prot.FCC") +
    theme_classic() +
    scale_fill_manual(values = c('tomato1','darkolivegreen2','powderblue', 'thistle')) +
    theme(legend.position = "none") +
    theme(
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25)
    ) +
    stat_summary(fun.data = "mean_sdl", geom = "pointrange", color = "black", position = position_dodge(0.75))
  
  plots[[paste0("base", position, "plot")]] <- base_plot
}

# Arrange the plots
arranged_plots <- ggarrange(
  plots$base2plot, plots$base3plot, plots$base4plot,
  plots$base5plot, plots$base6plot, plots$base7plot,
  plots$base8plot, plots$base9plot, plots$base10plot,
  plots$base11plot, ncol = 2, nrow = 5
)

# Save and print the arranged plots
pdf('/Users/evehallett/Documents/dissertation/to_submit/Figure5/all_positions_violins.pdf', width = 40, height = 60)
print(arranged_plots)
dev.off()
