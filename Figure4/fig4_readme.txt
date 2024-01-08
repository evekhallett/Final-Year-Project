The scripts do as follows:



1. get_acc_SE.py 

EC_SE_proteins.csv originally had 3 columns, denoting E. coli coding sequence accession numbers and their corresponding proteins 
and a column for orthologous protein names in Salmonella enterica (this was the output from BioCyc)

This script finds the S. enterica accession numbers from the corresponding 
protein names in the orthologs fasta and makes a column called 'SE_protein' 


2. get_orthologs_SE_EC.py 

SE.fa = unfiltered Salmonella enterica coding sequences from BioCyc

Uses the protein name/accession numbers from the S. enterica orthologs and extracts the sequences for these out of SE.fasta


3. make_ortholog_files.py

This script makes a fasta file for each orthologous sequence pair in E. coli and S. enterica. 

I then put these in a folder called 'unaligned_sequences'


4. MAFFT.py

Goes through each fasta in folder and aligns them, outputting fasta files in 'aligned_sequences' directory


5. get_conservation_gc3_EC_SE.py

Filters the aligned sequences (no internal stops, starts with a start codon, exclusively atcg, no indels which are not a multiple of 3)

Finds the conservation of GC3 and the absolute conservation at the third nucleotide sites


6. gc_conservation.R 

Plots the data
