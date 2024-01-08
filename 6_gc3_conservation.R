library(ggplot2)
library(ggpmisc)
setwd('/Users/evehallett/Documents/dissertation/to_submit/Figure4')

#ABSOLUTE CONSERVATION
general_conservation <- read.csv('/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/EC_salmonella_conservation_absolute.csv')
general_conservation_df = data.frame(general_conservation)
general_conservation_df$Conservation <- "Absolute"

#GC3 CONSERVATION
gc3 <- read.csv('/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/EC_salmonella_conservation_gc3.csv')
gc3df = data.frame(gc3)
gc3df$Conservation <- "GC3"

#MAKE ONE DF
combined_df <- rbind(gc3df, general_conservation_df)
combined_df$percentage_conservation <- round(combined_df$percentage_conservation, digits = 3)

#ADJUST COLOURS AND TITLES
my_colours <- c("GC3" = "deepskyblue3", "Absolute" = "brown3")
yname <- expression(paste("% Conservation between ", italic("E.coli "), " and ", italic("S. enterica")))

#PLOT
pdf("Fig4_0701.pdf")
plot3 <- ggplot(data = combined_df, aes(x = codon_number, y = percentage_conservation, colour = Conservation)) +
  stat_poly_line(formula = y ~ poly(x, 3, raw = TRUE), se = FALSE) +
  stat_poly_eq(
    formula = y ~ poly(x, 3, raw = TRUE),
    use_label(c("R2", "p")), 
    position = "identity",
    parse = TRUE
  ) +
  geom_point() +
  theme_classic() + 
  coord_cartesian(ylim = c(47, 100)) +
  scale_color_manual(values = my_colours, labels = c("Absolute", expression(GC[3]))) +
  ylab(yname) + xlab("Codon Position")

print(plot3)
dev.off()

