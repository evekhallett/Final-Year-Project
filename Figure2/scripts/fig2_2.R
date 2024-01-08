setwd("/Users/evehallett/Documents/dissertation/to_submit/Figure2/files/")
library(boot)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggrepel)
library(viridis)
perform_analysis <- function(file_path, col_name_prefix) {
  gc_df <- read.table(file_path, sep = ",", header = 1)
  
  bootstrap_results <- data.frame()
  boot.mean <- function(gc_df, i) {
    boot.mean <- mean(gc_df[i])
  }
  
  for (col_name in colnames(gc_df)) {
    result <- boot.mean(gc_df[[col_name]])
    boot_sd <- sd(boot(gc_df[[col_name]], boot.mean, R = 1000)$t)
    codon_number <- as.numeric(sub(col_name_prefix, "", col_name))
    bootstrap_results <- rbind(bootstrap_results, data.frame(Codon_number = codon_number, Mean = result[1], sd = boot_sd))
  }
  
  return(list(bootstrap_results = bootstrap_results))
}

#APPLY FUNCTION TO EACH CSV FILE, ADD RESULTS TO DATAFRAME AND PLOT BOOTSTRAPPED MEANS
gc1_analysis <- perform_analysis("ecoli_gc1_40.csv", "codon_")
gc1_results <- data.frame(gc = "gc1", gc1_analysis$bootstrap_results)
#gc1_plot <- gc1_analysis$plot+ labs(x= "Codon Number",y=("Mean GC1")) + ylim(0.3, 0.61)


gc2_analysis <- perform_analysis("ecoli_gc2_40.csv", "codon_")
gc2_results <- data.frame(gc = "gc2", gc2_analysis$bootstrap_results)
#gc2_plot <- gc2_analysis$plot+ labs(x= "Codon Number",y=("Mean GC2")) + ylim(0.3, 0.61)


gc3_analysis <- perform_analysis("ecoli_gc3_40.csv", "codon_")
gc3_results <- data.frame(gc = "gc3", gc3_analysis$bootstrap_results)
#gc3_plot <- gc3_analysis$plot+ labs(x= "Codon Number",y=("Mean GC3")) + ylim(0.3, 0.61)

all_results <- rbind(gc1_results, gc2_results, gc3_results)
print(all_results)



#MAKE DATAFRAME OF JUST THE BOOSTRAPPED MEANS
all_means <- unstack(all_results, Mean ~ gc, header=TRUE)
all_means$codon_number<-c(2:40)
print(all_means)


# USED THIS ONE
create_scatter_plot <- function(data, x_col, y_col, label_col, xlab, ylab) {
  pca <- prcomp(~get(x_col) + get(y_col), data)
  slp <- with(pca, rotation[2, 1] / rotation[1, 1])
  int <- with(pca, center[2] - slp * center[1])
  plot <- ggplot(data, aes_string(x = x_col, y = y_col, label = label_col, color = "codon_number")) +
    geom_point(shape = 19, size = 2) +
    labs(x = xlab, y = ylab, color="Codon Position") +
    geom_text_repel(hjust = -1, vjust = 1, size = 4, color = 'black') +
    theme_classic() +
    geom_abline(slope=slp, intercept =int, col ="brown")+
    stat_cor(label.y = 0.3, label.x = 1, method = "spearman") +
    annotate("text", x = 0.65, y = 0.28, label = bquote(italic("R")^2 ~ "=" ~ .(round(summary(lm(data[[y_col]] ~ data[[x_col]]))$r.squared, 3))), hjust = 1, vjust = -2) +
    annotate("text", x = 0.65, y = 0.25, label = bquote(italic("p") ~ "=" ~ .(formatC(summary(lm(data[[y_col]] ~ data[[x_col]]))$coef[2, 4], digits = 3))), hjust = 1, vjust = -2.5) +
    xlim(0.3, 0.65) + ylim(0.25, 0.65) +
    scale_color_viridis_c(option = "viridis") +
    scale_fill_viridis()
  return(plot)
}

# MAKE SCATTER PLOTS COMPARING GC1, 2 AND 3
plot2vs3 <- create_scatter_plot(all_means, "gc3", "gc2", "codon_number", expression("Mean GC"[3]), expression("Mean GC"[2]))
plot1vs2 <- create_scatter_plot(all_means, "gc2", "gc1", "codon_number", expression("Mean GC"[2]), expression("Mean GC"[1]))
plot3vs1 <- create_scatter_plot(all_means, "gc3", "gc1", "codon_number", expression("Mean GC"[3]), expression("Mean GC"[1]))

pdf('Figure2_0701.pdf', width = 12, height = 10)
plots <- ggarrange(plot2vs3, plot1vs2, plot3vs1, 
          labels = c("a", "b", "c"),
          ncol = 2, nrow = 2)
print(plots)
dev.off()


