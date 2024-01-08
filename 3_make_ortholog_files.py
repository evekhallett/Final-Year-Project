import pandas as pd
from Bio import SeqIO
import re
import os

os.chdir("/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/unaligned_sequences")
id_file = pd.read_csv("/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/EC_SA_proteins.csv")
id_file = id_file.dropna()
SE_fasta = '/Users/evehallett/Documents/dissertation/salmonella_ecoli_orthologs/SE_ortholog_seqs.fa'
EC_fasta = '/Users/evehallett/Documents/dissertation/salmonella_ecoli_orthologs/EC_ortholog_seqs.fa'

SE_id = id_file['SE_accession'].astype(str).str.strip().tolist()
EC_id = id_file['EC_accession'].astype(str).str.strip().tolist()

ids = dict(zip(SE_id, EC_id))
def is_valid_sequence(record):
    # CHECK SEQUENCE MULTIPLE OF 3 BASES
    if len(record.seq) % 3 != 0:
        return False

    # CHECK SEQUENCE ONLY CONTAINS ATCG
    if set(record.seq.upper()) - set("ATCG"):
        return False

    # CHECK SEQUENCE STARTS WITH "NTG" AND ENDS IN A STOP
    start_motif = re.match(r'[A-Za-z][Tt][Gg]', str(record.seq))
    end_motif = re.search(r'(TAA|TGA|TAG)[^ATCG]*$', str(record.seq))

    # EXCLUDE IF HAS INTERNAL STOP
    for i in range(4,len(record.seq)-2):
        codon = str(record.seq[i:i+2])
        if re.search(r'(TAA|TGA|TAG)', codon):
            print(f"Internal stop codon {codon} found at position {i} for {record.id}")
            return False

    if not start_motif or not end_motif:
        return False

    return True


def extract_sequences_by_id(fasta_file, locus_id):
    matching_sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        match = re.search(r'\b{}\b'.format(re.escape(locus_id)), record.description)
        if match and is_valid_sequence(record):
            matching_sequences.append(record)
    return matching_sequences

def extract_sequences(SE_fasta, EC_fasta, ids):
    for locus_id, output_file in ids.items():
        #NAME FASTA AFTER S.ENTERICA ACCESSION
        output_file = '{}.fasta'.format(locus_id)
        SE_sequences = extract_sequences_by_id(SE_fasta, locus_id)
        EC_sequences = extract_sequences_by_id(EC_fasta, ids[locus_id])

        print(f"Number of SE sequences for {locus_id}: {len(SE_sequences)}")
        print(f"Number of EC sequences for {locus_id}: {len(EC_sequences)}")

        # Check the validity of sequences before combining
        valid_SE_sequences = [seq for seq in SE_sequences if is_valid_sequence(seq)]
        valid_EC_sequences = [seq for seq in EC_sequences if is_valid_sequence(seq)]

        combined_sequences = valid_EC_sequences + valid_SE_sequences

        if combined_sequences:
            with open(output_file, 'w') as out_fasta:
                SeqIO.write(combined_sequences, out_fasta, 'fasta')

extract_sequences(SE_fasta, EC_fasta, ids)
