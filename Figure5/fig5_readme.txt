The scripts do as follows:

1. goodman_nt3.py

Extracts the Prot.FCC values and third nucleotides for codons 2-11 from coding sequences, input file is data from Goodman, Church and Kosuri (2013)
Output is goodman_nt3.csv, which is also used for figure 6.

2. goodman_violins.R

Plots the frequency of each base (with median and standard deviation) at the third nucleotide site for each codon position in the coding sequences extracted.


3. anova+tukey.R

Runs anova analysis between the Prot.FCC values and the nucleotides present at the third site of codons 2-11. Extracts and plots F statistic vs codon postion. 

Points in the anova plot are meant to be colour coded based on the p value, however they were all below 0.001, hence in the output they are the same colour.

Also runs Tukey post-hoc tests between each base combination with Prot.FCC. Results are then plotted separately based on codon position (Figure S3) 
