library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggpubr)
setwd('/Users/evehallett/Documents/dissertation/to_submit/Figure6/files/')
all_codons<- read.csv('goodman_nt3.csv')

base_columns <- c("Position_2", "Position_3", "Position_4", "Position_5", "Position_6", "Position_7", "Position_8", "Position_9", "Position_10", "Position_11")
#WANT TO FIND THE MEAN FREQUENCY OF EACH NUCLEOTIDE FOR CODON POSITIONS 2-10 AND 11
means_df <- data.frame(Name = character(), 
                       Mean_A_2_to_10 = numeric(), Mean_A_11 = numeric(), 
                       Mean_T_2_to_10 = numeric(), Mean_T_11 = numeric(), 
                       Mean_G_2_to_10 = numeric(), Mean_G_11 = numeric(), 
                       Mean_C_2_to_10 = numeric(), Mean_C_11 = numeric(), 
                       stringsAsFactors = FALSE)

for (i in 1:nrow(all_codons)) {
  name <- all_codons$Name[i]
  
  #A
  mean_A_2_to_10 <- mean(all_codons[i, base_columns[1:9]] == 'A')
  mean_A_11 <- mean(all_codons[i, "Position_11"] == 'A')
  
  #T
  mean_T_2_to_10 <- mean(all_codons[i, base_columns[1:9]] == 'T')
  mean_T_11 <- mean(all_codons[i, "Position_11"] == 'T')
  
  #G
  mean_G_2_to_10 <- mean(all_codons[i, base_columns[1:9]] == 'G')
  mean_G_11 <- mean(all_codons[i, "Position_11"] == 'G')
  
  #C
  mean_C_2_to_10 <- mean(all_codons[i, base_columns[1:9]] == 'C')
  mean_C_11 <- mean(all_codons[i, "Position_11"] == 'C')
  
  #APPEND MEANS TO SINGLE DF
  means_df <- rbind(means_df, data.frame(Name = name, 
                                         Mean_A_2_to_10 = mean_A_2_to_10, Mean_A_11 = mean_A_11, 
                                         Mean_T_2_to_10 = mean_T_2_to_10, Mean_T_11 = mean_T_11,
                                         Mean_G_2_to_10 = mean_G_2_to_10, Mean_G_11 = mean_G_11,
                                         Mean_C_2_to_10 = mean_C_2_to_10, Mean_C_11 = mean_C_11))
}

print(means_df)

#PUT MEANS DF INTO CSV
write.csv(means_df, file = 'means.csv', row.names = FALSE)

#WILCOXON TEST
perform_wilcoxon_test <- function(group1, base1, group2, base2) {
  result <- wilcox.test(group1, group2)
  return(data.frame(Test_Type = paste("Mean_", base1, "_2_to_10 vs Mean_", base2, "_11", sep = ""), 
                    W_Statistic = result$statistic, P_Value = result$p.value))
}
#WILCOXON TEST RESULTS DF
wilcoxon_results_df <- data.frame(Test_Type = character(), W_Statistic = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

#TESTS BETWEEN EACH BASE COMBINATION
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_A_2_to_10, 'A', means_df$Mean_A_11, "A"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_A_2_to_10, 'A', means_df$Mean_T_11, "T"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_A_2_to_10, 'A', means_df$Mean_C_11, "C"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_A_2_to_10, 'A', means_df$Mean_G_11, "G"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_T_2_to_10, 'T', means_df$Mean_A_11, "A"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_T_2_to_10, 'T', means_df$Mean_T_11, "T"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_T_2_to_10, 'T', means_df$Mean_C_11, "C"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_T_2_to_10, 'T', means_df$Mean_G_11, "G"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_C_2_to_10, 'C', means_df$Mean_A_11, "A"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_C_2_to_10, 'C', means_df$Mean_T_11, "T"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_C_2_to_10, 'C', means_df$Mean_C_11, "C"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_C_2_to_10, 'C', means_df$Mean_G_11, "G"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_G_2_to_10, 'G', means_df$Mean_A_11, "A"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_G_2_to_10, 'G', means_df$Mean_T_11, "T"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_G_2_to_10, 'G', means_df$Mean_C_11, "C"))
wilcoxon_results_df <- rbind(wilcoxon_results_df, perform_wilcoxon_test(means_df$Mean_G_2_to_10, 'G', means_df$Mean_G_11, "G"))

print(wilcoxon_results_df)


A_df <- wilcoxon_results_df[grep("Mean_A_2_to_10", wilcoxon_results_df$Test_Type), ]
A_df$Comparison1 <- c("A","T","C","G")
A_df$Comparison1 <- factor(A_df$Comparison1)
A_df$Comparison2 <- c("A","A","A","A")
A_df$Comparison2 <- factor(A_df$Comparison2)

T_df <- wilcoxon_results_df[grep("Mean_T_2_to_10", wilcoxon_results_df$Test_Type), ]
print(T_df)
T_df$Comparison1 <- c("A","T","C","G")
T_df$Comparison1 <- factor(T_df$Comparison1)
T_df$Comparison2 <- c("T","T","T","T")
T_df$Comparison2 <- factor(T_df$Comparison2)

G_df <- wilcoxon_results_df[grep("Mean_G_2_to_10", wilcoxon_results_df$Test_Type), ]
G_df$Comparison1 <- c("A","T","C","G")
G_df$Comparison1 <- factor(G_df$Comparison1)
G_df$Comparison2 <- c("G","G","G","G")
G_df$Comparison2 <- factor(G_df$Comparison2)

C_df <- wilcoxon_results_df[grep("Mean_C_2_to_10", wilcoxon_results_df$Test_Type), ]
C_df$Comparison1 <- c("A","T","C","G")
C_df$Comparison1 <- factor(C_df$Comparison1)
C_df$Comparison2 <- c("C","C","C","C")
C_df$Comparison2 <- factor(C_df$Comparison2)


all_codons_df <- rbind(A_df,T_df,C_df,G_df)
print(all_codons_df)
desired_order <- c('A', 'T', 'C', 'G')

# AUTOMATICALLY PUTS BASES IN ALPHABETICAL ORDER SO USE THIS TO CUSTOMISE
all_codons_df$Comparison1 <- factor(all_codons_df$Comparison1, levels = desired_order)

pdf('wilcoxon_plot.pdf')
#LINE PLOT WHICH GROUPS EACH SET OF COMPARISONS BY THE NUCLEOTIDE AT THE 3RD POSITION CODONS 2-10
plot_wilcoxon <- ggplot(all_codons_df, aes(x = Comparison1, y = W_Statistic, colour = Comparison2, group = Comparison2)) +
  geom_point() + 
  geom_line()+
  scale_x_discrete(name = 'Nucleotide 3, Codon Position 11') +
  scale_y_continuous(name = 'W Statistic') +
  theme_classic() + labs(color = "Nucleotide 3, Codon Positions 2-10")

print(plot_wilcoxon)
dev.off()
