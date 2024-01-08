1. Calculating the GC1, 2 and 3 content of the 5' ends of CDS

GC1, GC2 and GC3.py calculates the GC1, 2 and 3 content of the 5'-end codons in the 
supplied fasta file containing E. coli CDS. 

Filters the CDS so ones which start with NTG, end in a stop codon, contain only ATCG and 
have no internal stop codons are included in the analysis.

The outputs are ecoli_gc1_40.csv, ecoli_gc3_40.csv and ecoli_gc3_40.csv which are vectors 
of GC content for the first 40 codons. 

2. Plotting the bootstrapped mean values - Fig1.R

Uses the boot package to bootstrapped the mean GC content for each of the csv files and 
plot them with standard error bars, with a separate plot for each nucleotide position.

3. Fig2.R 

PCAs of the bootstrapped mean GC1, 2 and 3 content