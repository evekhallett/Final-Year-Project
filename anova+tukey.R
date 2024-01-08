library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ggpmisc)
setwd('/Users/evehallett/Documents/dissertation/to_submit/Figure5/')
df<- read.csv('/Users/evehallett/Documents/dissertation/to_submit/Figure5/goodman_nt3.csv')
statistics_df <- data.frame(Position = character(), F_statistic = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

#ANOVA TESTS
oneway2 <- aov(Prot.FCC ~ Position_2, data = df)
oneway3 <- aov(Prot.FCC ~ Position_3, data = df)
oneway4 <- aov(Prot.FCC ~ Position_4, data = df)
oneway5 <- aov(Prot.FCC ~ Position_5, data = df)
oneway6 <- aov(Prot.FCC ~ Position_6, data = df)
oneway7 <- aov(Prot.FCC ~ Position_7, data = df)
oneway8 <- aov(Prot.FCC ~ Position_8, data = df)
oneway9 <- aov(Prot.FCC ~ Position_9, data = df)
oneway10 <- aov(Prot.FCC ~ Position_10, data = df)
oneway11 <- aov(Prot.FCC ~ Position_11, data = df)

summary_oneway2 <- summary(oneway2)
summary_oneway3 <- summary(oneway3)
summary_oneway4 <- summary(oneway4)
summary_oneway5 <- summary(oneway5)
summary_oneway6 <- summary(oneway6)
summary_oneway7 <- summary(oneway7)
summary_oneway8 <- summary(oneway8)
summary_oneway9 <- summary(oneway9)
summary_oneway10 <- summary(oneway10)
summary_oneway11 <- summary(oneway11)

#ADDS EACH TEST RESULT TO A SINGLE DATAFRAME
statistics_df <- rbind(statistics_df, data.frame(Position = "2", F_statistic = summary_oneway2[[1]]$`F value`, p_value = summary_oneway2[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "3", F_statistic = summary_oneway3[[1]]$`F value`, p_value = summary_oneway3[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "4", F_statistic = summary_oneway4[[1]]$`F value`, p_value = summary_oneway4[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "5", F_statistic = summary_oneway5[[1]]$`F value`, p_value = summary_oneway5[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "6", F_statistic = summary_oneway6[[1]]$`F value`, p_value = summary_oneway6[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "7", F_statistic = summary_oneway7[[1]]$`F value`, p_value = summary_oneway7[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "8", F_statistic = summary_oneway8[[1]]$`F value`, p_value = summary_oneway8[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "9", F_statistic = summary_oneway9[[1]]$`F value`, p_value = summary_oneway9[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "10", F_statistic = summary_oneway10[[1]]$`F value`, p_value = summary_oneway10[[1]]$"Pr(>F)"[1]))
statistics_df <- rbind(statistics_df, data.frame(Position = "11", F_statistic = summary_oneway11[[1]]$`F value`, p_value = summary_oneway11[[1]]$"Pr(>F)"[1]))
statistics_df <- na.omit(statistics_df)

#MAKES SURE CODON POSITIONS ARE PLOTTED CHRONOLOGICALLY
statistics_df$Position <- as.integer(statistics_df$Position)

#ANOVA PLOT
pdf('W-stat-plot.pdf')
statistics_df$Significance <- ifelse(statistics_df$p_value < 0.001, "p < 0.001",
                                     ifelse(statistics_df$p_value < 0.01, "p < 0.01", "p < 0.05"))


plot <- ggplot(data = statistics_df, aes(x = Position, y = F_statistic, color = Significance)) +
  geom_point(shape = 18) +
  stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +
  stat_poly_eq(
    formula = y ~ poly(x, 2, raw = TRUE),
    use_label(c("R2", "p")),
    position = "identity",
    parse = TRUE
  ) +
  scale_color_manual(values = c("p < 0.001" = "blue", "p < 0.01" = "green", "p < 0.05" = "orange")) +
  theme_classic() +
  labs(x = 'Codon Position', y = 'F Statistic')

print(plot)
dev.off()

#TUKEY POST HOC
tukey2 <- TukeyHSD(oneway2)
tukey3 <- TukeyHSD(oneway3)
tukey4 <- TukeyHSD(oneway4)
tukey5 <- TukeyHSD(oneway5)
tukey6 <- TukeyHSD(oneway6)
tukey7 <- TukeyHSD(oneway7)
tukey8 <- TukeyHSD(oneway8)
tukey9 <- TukeyHSD(oneway9)
tukey10 <- TukeyHSD(oneway10)
tukey11 <- TukeyHSD(oneway11)

#TUKEY PLOTS
plot(tukey2, las=1, col = 'brown')
plot(tukey3, las=1, col = 'brown')
plot(tukey4, las=1, col = 'brown')
plot(tukey5, las=1, col = 'brown')
plot(tukey6, las=1, col = 'brown')
plot(tukey7, las=1, col = 'brown')
plot(tukey8, las=1, col = 'brown')
plot(tukey9, las=1, col = 'brown')
plot(tukey10, las=1, col = 'brown')
plot(tukey11, las=1, col = 'brown')

