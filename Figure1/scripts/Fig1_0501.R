setwd("/Users/evehallett/Documents/dissertation/to_submit/Figure1/files/")
library(boot)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(gridExtra)
library(robustbase)
library(ggtext)

#FINDS THE POINTS OF INFLECTION 
find_inflection_point <- function(x, y) {
  derivative <- diff(y) / diff(x)
  
  #WHEN THE DIFFERENCE IN THE DERIVATIVES IS GREATER THAN 0, I.E. CHANGES SIGN
  sign_change_indices <- which(diff(derivative > 0) != 0)
  
  # RETURN 1ST POINT OF INFLECTION
  if (length(sign_change_indices) > 0) {
    return(c(x[sign_change_indices[1]], y[sign_change_indices[1]]))
  } else {
    return(NULL)
  }
}


perform_analysis <- function(file_path, col_name_prefix, polynomial_degree) {
  # IMPORT DATAFRAME
  gc_df <- data.frame()
  gc_table <- read.table(file_path, sep = ",", header = 1)
  gc_df <- data.frame(gc_table)
  
  bootstrap_results <- data.frame()
  # BOOTSTRAPPING FUNCTION
  boot.mean <- function(gc_df, i) {
    boot.mean <- mean(gc_df[i])
  }
  
  for (col_name in colnames(gc_df)) {
    result <- boot.mean(gc_df[[col_name]])
    boot_sd <- sd(boot(gc_df[[col_name]], boot.mean, R = 1000)$t)
    codon_number <- as.numeric(sub(col_name_prefix, "", col_name))
    bootstrap_results <- rbind(bootstrap_results, data.frame(Codon_number = codon_number, Mean = result[1], sd = boot_sd))
  }
  
  #SCATTER WITH ERROR BARS
  plot_data <- ggplot_build(ggplot(bootstrap_results, aes(x = Codon_number, y = Mean)) +
                              geom_point(col = "grey") +
                              geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), col = "grey") +
                              labs(x = "Codon Number", y = "Mean GC") +
                              theme_classic())$data[[1]]
  
  # POLYNOMIAL FIT
  polynomial_fit <- lm(Mean ~ poly(Codon_number, degree = polynomial_degree), data = bootstrap_results)
  polynomial_pred <- predict(polynomial_fit, newdata = data.frame(Codon_number = plot_data$x))
  ssr_polynomial <- sum((polynomial_pred - plot_data$y)^2)
  sst_polynomial <- sum((plot_data$y - mean(plot_data$y))^2)
  # R^2 AND P VALUE OF POLYNOMIAL FIT
  r_squared_polynomial <- 1 - (ssr_polynomial / sst_polynomial)
  p_value_polynomial <- summary(polynomial_fit)$coefficients[, "Pr(>|t|)"]
  
  #APPLY INFLECTION POINT FUNCTION
  inflection_point <- find_inflection_point(plot_data$x, polynomial_pred)
  
  #PRINT CODON POSITION WHERE GRADIENT SWITCHES FROM POSITIVE TO NEGATIVE
  if (!is.null(inflection_point)) {
    cat("Point where gradient switches from positive to negative:", inflection_point[1], "\n")
    
    #PLOTS POLYNOMIAL WITH R^2, P VALUE AND POINT OF INFLECTION
    plot <- ggplot() + 
      geom_point(data = bootstrap_results, aes(x = Codon_number, y = Mean), col = "grey") +
      geom_errorbar(data = bootstrap_results, aes(x = Codon_number, ymin = Mean - sd, ymax = Mean + sd), col = "grey") +
      geom_line(data = data.frame(Codon_number = plot_data$x, Mean = polynomial_pred),
                aes(x = Codon_number, y = Mean), col = "deepskyblue4") +
      annotate("point", x = inflection_point[1], y = inflection_point[2], col = "brown", size = 3) +
      labs(x = "Codon Position", y = "Mean GC") +
      theme_classic()+
      annotate("text", x = 40, y = 0.325, label = bquote(italic("R")^2 ~ "=" ~ .(round(r_squared_polynomial, 3))), hjust = 1, vjust = -2) +
      annotate("text", x = 40, y = 0.30, label = bquote(italic("p") ~ "=" ~ .(formatC(p_value_polynomial, digits = 3))), hjust = 1, vjust = -2.5)
    
  } else {
    cat("No point where gradient switches from positive to negative.\n")
    #PLOT WITHOUT POINT OF INFLECTION
    plot <- ggplot() + 
      geom_point(data = bootstrap_results, aes(x = Codon_number, y = Mean), col = "grey") +
      geom_errorbar(data = bootstrap_results, aes(x = Codon_number, ymin = Mean - sd, ymax = Mean + sd), col = "grey") +
      geom_line(data = data.frame(Codon_number = plot_data$x, Mean = polynomial_pred),
                aes(x = Codon_number, y = Mean), col = "deepskyblue4") +
      labs(x = "Codon Position", y = "Mean GC") +
      theme_classic()+
      annotate("text", x = 40, y = 0.315, label = bquote(italic("R")^2 ~ "=" ~ .(round(r_squared_polynomial, 3))), hjust = 1, vjust = -2) +
      annotate("text", x = 40, y = 0.30, label = bquote(italic("p") ~ "=" ~ .(formatC(p_value_polynomial, digits = 3))), hjust = 1, vjust = -2.5)
    
  }
  
  return(list(plot = plot, bootstrap_results = bootstrap_results))
}


#CHANGE TO SEE WHAT FITS BEST
polynomial_degree <- 4

#GC1
gc1_analysis <- perform_analysis("ecoli_gc1_40.csv", "codon_", polynomial_degree)
gc1_results <- data.frame(gc = "gc1", gc1_analysis$bootstrap_results)
gc1_plot <- gc1_analysis$plot + labs(x= "Codon Position",y = expression("Mean GC"[1])) + ylim(0.3, 0.62)
print(gc1_plot)

#GC2
gc2_analysis <- perform_analysis("ecoli_gc2_40.csv", "codon_", polynomial_degree)
gc2_results <- data.frame(gc = "gc2", gc2_analysis$bootstrap_results)
gc2_plot <- gc2_analysis$plot + labs(x= "Codon Position",y = expression("Mean GC"[2])) + ylim(0.3, 0.62)

#GC3
gc3_analysis <- perform_analysis("ecoli_gc3_40.csv", "codon_", polynomial_degree)
gc3_results <- data.frame(gc = "gc3", gc3_analysis$bootstrap_results)
gc3_plot <- gc3_analysis$plot + labs(x= "Codon Position", y = expression("Mean GC"[3])) + ylim(0.3, 0.62)

#PLOT TOGETHER
pdf('fig1_0701.pdf')
plots <- ggarrange(gc1_plot, gc2_plot, gc3_plot, 
                   labels = c("a", "b", "c"),
                   ncol = 2, nrow = 2)
print(plots)
dev.off()