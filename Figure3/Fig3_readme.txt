1. translation.py

Extracts the 2nd to the 40th codons in the CDS of the supplied E. coli fasta file.

For each amino acid, the frequency and relative frequency of each codon is calculated for 
each codon position. These are outputted in a unique csv for each amino acid, deposited 
in the 'amino_acids_final' directory.

2. all_amino_acids.R

Creates line plots of the percentage relative codon usage frequencies, with one plot per
amino acid. Colour coding based on GC content (general and GC3).