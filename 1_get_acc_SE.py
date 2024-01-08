from Bio import SeqIO
import pandas as pd

#CSV FILE ORIGINALLY HAS COLUMNS: "EC_protein". "EC_accession". "SE_protein", script ultimately adds "SE_accession"
csv_file = "/Users/evehallett/Documents/dissertation/to_submit/Figure3/files/EC_SA_proteins.csv"
#FASTA WITH
fasta_file = "/Users/evehallett/Documents/dissertation/to_submit/Figure3/files/SE.fasta"

df = pd.read_csv(csv_file)
protein_names = df['SE_protein'].tolist()

#DICTIONARY TO APPEND S. ENTERICA ACCESSION NUMBERS TO PROTEIN NAMES
accession_mapping = {}

with open(fasta_file, 'r'):
    for record in SeqIO.parse(fasta_file, 'fasta'):
        header = record.description
        #EXTRACTS SECOND ELEMENT BETWEEN PIPES
        accession = header.split('|')[1] if '|' in header else ''
        #protein name is the last part of the accession after the last | - make sure this is first as will search for string inside protein names
        protein_name = header.split('|')[-1].strip()
        #IF PROTEIN IN LIST, WILL ADD ACCESSION TO THE DICT
        if protein_name in protein_names:
            accession_mapping[protein_name] = accession

#ADDS THE ACCESSION NUMBER CORRECT ROW IN NEW COLUMN FOR S. ENTERICA ACCESSION NUMBERS
df['SE_accession'] = df['SE_protein'].map(accession_mapping)
df.to_csv(csv_file, index=False)
