library(ggplot2)
library(ggpubr)

folder_path <- "/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/"
setwd(folder_path)
#CSV FILES NAMED AFTER AMINO ACID, MAKES LIST 
file_list <- list.files(folder_path, pattern = "\\.csv", full.names = TRUE)


for (file_name in file_list) {
  df <- read.csv(file_name)
  df <- data.frame(df)
  #FIRST 20 CODONS
  df <- df[df$codon_number <= 20, ]
  #NAME DATAFRAME AFTER AMINO ACID
  df_name <- gsub("\\.csv", "", basename(file_name))
  assign(df_name, df)
}


alanine$codon <- as.factor(alanine$codon)
alanine_plot <- ggplot(alanine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("GCA" = 'red', "GCC" = 'darkorange1', "GCG"='darkolivegreen3', "GCT"= 'deepskyblue3')) +
  labs(color = "Alanine Codon", linetype = "Alanine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("GCA" = 'longdash', "GCC" = 'solid', "GCG"='solid', "GCT"= 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

arginine$codon <- as.factor(arginine$codon)
arginine_plot <- ggplot(arginine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("AGA" = "deeppink3", "AGG" = "darkgreen", "CGA" = 'red', "CGC" = 'darkorange1', "CGG"='darkolivegreen3', "CGT"= 'deepskyblue3')) +
  labs(color = "Arginine Codon", linetype = "Arginine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("AGA" = "twodash", "AGG" = "longdash", "CGA" = 'longdash', "CGC" = 'solid', "CGG"='solid', "CGT"= 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

asparagine$codon <- as.factor(asparagine$codon)
asparagine_plot <- ggplot(asparagine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("AAC" = 'darkorange1', "AAT" = 'deepskyblue3')) +
  labs(color = "Asparagine Codon", linetype = "Asparagine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("AAC" = 'twodash', "AAT" = 'dotted'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

aspartate$codon <- as.factor(aspartate$codon)
aspartate_plot <- ggplot(aspartate, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("GAT" = 'deepskyblue3', "GAC" = 'darkorange1')) +
  labs(color = "Aspartate Codon", linetype = "Aspartate Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("GAC"='longdash', "GAT"= 'twodash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

cystine$codon <- as.factor(cystine$codon)
cystine_plot <- ggplot(cystine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("TGT" = 'deepskyblue3', "TGC" = 'darkorange1')) +
  labs(color = "Cystine Codon", linetype = "Cystine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("TGT" = 'twodash', "TGC" = 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

glutamate$codon <- as.factor(glutamate$codon)
glutamate_plot <- ggplot(glutamate, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("GAA" = 'red', "GAG" = 'darkolivegreen3')) +
  labs(color = "Glutamate Codon", linetype = "Glutamate Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("GAA"='twodash', "GAG"= 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

glutamine$codon <- as.factor(glutamine$codon)
glutamine_plot <- ggplot(glutamine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("CAA" = 'red', "CAG" = 'darkolivegreen3')) +
  labs(color = "Glutamine Codon", linetype = "Glutamine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("CAA" = 'twodash', "CAG" = 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

glycine$codon <- as.factor(glycine$codon)
glycine_plot <- ggplot(glycine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("GGA" = 'red', "GGC" = 'darkorange1', "GGG"='darkolivegreen3', "GGT"= 'deepskyblue3')) +
  labs(color = "Glycine Codon", linetype = "Glycine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("GGA" = 'longdash', "GGC" = 'solid', "GGG"='solid', "GGT"= 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

histidine$codon <- as.factor(histidine$codon)
histidine_plot <- ggplot(histidine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("CAT" = 'deepskyblue3', "CAC" = 'darkorange1')) +
  labs(color = "Histidine Codon", linetype = "Histidine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("CAT" = 'twodash', "CAC" = 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

isoleucine$codon <- as.factor(isoleucine$codon)
isoleucine_plot <- ggplot(isoleucine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("ATT" = 'deepskyblue3', "ATC" = 'darkorange1', 'ATA' = 'red')) +
  labs(color = "Isoleucine Codon", linetype = "Isoleucine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("ATT" = 'dotted', "ATC" = 'twodash', 'ATA' = 'dotted'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

leucine$codon <- as.factor(leucine$codon)
leucine_plot <- ggplot(leucine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("TTA"='deeppink3', "TTG"= 'darkgreen',"CTA" = "red", "CTC" = "darkorange1", "CTG" = 'darkolivegreen3', "CTT" = 'deepskyblue3')) +
  labs(color = "Leucine Codon", linetype = "Leucine Codon") + 
  xlab("Codon Positon") + ylab("Relative Frequency") +
  theme_classic()+
  scale_linetype_manual(values = c("TTA"='dotted', "TTG"= 'twodash',"CTA" = "twodash", "CTC" = "longdash", "CTG" = 'longdash', "CTT" = 'twodash'))+
  ylim(0,1)+ theme(aspect.ratio=1)

lysine$codon <- as.factor(lysine$codon)
lysine_plot <- ggplot(lysine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("AAA" = 'red', "AAG" = 'darkolivegreen3')) +
  labs(color = "Lysine Codon", linetype = "Lysine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("AAA" = 'dotted', "AAG" = 'twodash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

phenylalanine$codon <- as.factor(phenylalanine$codon)
phenylalanine_plot <- ggplot(phenylalanine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("TTT" = 'deepskyblue3', "TTC" = 'darkorange1')) +
  labs(color = "Phenylalanine Codon", linetype = "Phenylalanine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("TTT" = 'dotted', "TTC" = 'twodash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

proline$codon <- as.factor(proline$codon)
proline_plot <- ggplot(proline, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("CCA" = 'red', "CCC" = 'darkorange1', "CCG"='darkolivegreen3', "CCT"= 'deepskyblue3')) +
  labs(color = "Proline Codon", linetype = "Proline Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("CCA" = 'longdash', "CCC" = 'solid', "CCG"='solid', "CCT"= 'longdash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

serine$codon <- as.factor(serine$codon)
serine_plot <- ggplot(serine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("AGC" = "brown", "AGT" = "blue", "TCA" = 'red', "TCC" = 'darkorange1', "TCG"='darkolivegreen3', "TCT"= 'deepskyblue3')) +
  labs(color = "Serine Codon", linetype = "Serine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("AGC" = "longdash", "AGT" = "twodash", "TCA" = 'twodash', "TCC" = 'longdash', "TCG"='longdash', "TCT"= 'twodash'))+
  ylim(0,1)+ theme(aspect.ratio=1)

threonine$codon <- as.factor(threonine$codon)
threonine_plot <- ggplot(threonine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("ACA" = 'red', "ACC" = 'darkorange1', "ACG"='darkolivegreen3', "ACT"= 'deepskyblue3')) +
  labs(color = "Threonine Codon", linetype = "Threonine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("ACA" = 'twodash', "ACC" = 'longdash', "ACG"='longdash', "ACT"= 'twodash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

tyrosine$codon <- as.factor(tyrosine$codon)
tyrosine_plot <- ggplot(tyrosine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("TAT" = 'deepskyblue3', "TAC" = 'darkorange1')) +
  labs(color = "Tyrosine Codon", linetype = "Tyrosine Codon") + 
  xlab("Codon Position") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("TAT" = 'dotted', "TAC" = 'twodash'))+ 
  ylim(0,1)+ theme(aspect.ratio=1)

valine$codon <- as.factor(valine$codon)
valine_plot <- ggplot(valine, aes(x = codon_number, y = relative_frequency, color = codon, linetype = codon)) +
  geom_line() +
  scale_color_manual(values = c("GTA" = 'red', "GTC" = 'darkorange1', "GTG"='darkolivegreen3', "GTT"= 'deepskyblue3')) +
  labs(color = "Valine Codon", linetype = "Valine Codon") + 
  xlab("Codon Number") + ylab("Relative Frequency") +
  theme_classic() +
  scale_linetype_manual(values = c("GTA" = 'twodash', "GTC" = 'longdash', "GTG"='longdash', "GTT"= 'twodash'))+ 
  ylim(0,1) + theme(aspect.ratio=1)

pdf("all_amino_acids.pdf", width = 20, height = 30)
plot <- ggarrange(arginine_plot, leucine_plot, serine_plot, alanine_plot, glycine_plot, proline_plot, threonine_plot, valine_plot, isoleucine_plot, asparagine_plot, aspartate_plot, cystine_plot, histidine_plot, phenylalanine_plot, tyrosine_plot, glutamate_plot, glutamine_plot, lysine_plot, 
          labels = c("a", "b","c", "d", "e", "f", "g","h","i","j","k","l","m","n","o","p","q", "r"),
          ncol = 3, nrow = 6)
print(plot)
dev.off()


pdf("amino_acids_no_6fold2.pdf", width = 15, height = 50)
plot_no6fold <- ggarrange(alanine_plot, glycine_plot, proline_plot, threonine_plot, valine_plot, isoleucine_plot, glutamate_plot, glutamine_plot, lysine_plot, asparagine_plot, aspartate_plot, cystine_plot, histidine_plot, phenylalanine_plot, tyrosine_plot, 
                  labels = c("a", "b","c", "d", "e", "f", "g","h","i","j","k","l","m","n","o"),
                  ncol = 2, nrow = 8)
print(plot_6fold)
dev.off()

pdf("6fold.pdf", width = 10, height = 10)
plot_no6fold <- ggarrange(arginine_plot, leucine_plot, serine_plot, 
                          labels = c("a", "b","c"),
                          ncol = 2, nrow = 2)
print(plot_6fold)
dev.off()

pdf("6fold_cystine.pdf", width = 10, height = 10)
plot_no6fold_cys <- ggarrange(arginine_plot, leucine_plot, serine_plot, cystine_plot,
                          labels = c("a", "b", "c", "d"),
                          ncol = 2, nrow = 2)
print(plot_no6fold_cys)
dev.off()
