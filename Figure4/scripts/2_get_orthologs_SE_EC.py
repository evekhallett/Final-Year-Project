from Bio import SeqIO
import pandas as pd

df = pd.read_csv("/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/EC_SA_proteins.csv")

fasta_file = "/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/SE.fasta"
output_file = "/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/SE_ortholog_seqs.fa"
protein_names = df['SE_protein'].tolist()

"""RUN SCRIPT AGAIN AND SWITCH TO THIS TO DO THE SAME FOR E. COLI
fasta_file = "/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/EC.fasta"
output_file = "/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/EC_ortholog_seqs.fa"
protein_names = df['EC_protein'].tolist()
"""

#USE ACCESSION NUMBERS TO PUT ALL S. ENTERICA OR ALL E. COLI ORTHOLOGS INTO A FASTA
with open(output_file, 'w') as out_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        header = record.description
        if any(f'|{str(name)}' in header for name in protein_names):
            SeqIO.write(record, out_file, 'fasta')
