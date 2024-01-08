from Bio import SeqIO
from collections import defaultdict, Counter
import pandas as pd
from Bio.SeqUtils import GC

def process_amino_acid(fasta_file, amino_acid_dict, output_filename):
    dicts_codons = defaultdict(Counter)

    for record in SeqIO.parse(fasta_file, "fasta"):
        coding_seq = str(record.seq)

        for n in range(5, 121, 3):
            cdn = coding_seq[n:n+3]
            pos = 1 + n//3
            dict_entry = f"{pos}_{cdn}"
            if cdn in amino_acid_dict:
                dicts_codons[dict_entry][cdn] += 1

    sorted_dict = dict(sorted(dicts_codons.items()))

    amino_acid_df = pd.DataFrame(sorted_dict.items(), columns=['codon_position','frequency'])
    amino_acid_df['frequency'] = amino_acid_df['frequency'].astype(str).replace('[^0-9]', '', regex=True)
    amino_acid_df[['codon_number', 'codon']] = amino_acid_df['codon_position'].str.split('_', 1, expand=True)
    amino_acid_df.drop('codon_position', axis=1, inplace=True)
    amino_acid_df['codon_number'] = pd.to_numeric(amino_acid_df['codon_number'])
    amino_acid_df['frequency'] = pd.to_numeric(amino_acid_df['frequency'])

    #COLUMN NAMES
    amino_acid_df = amino_acid_df[['codon_number', 'codon', 'frequency']]

    #SORT EVERYTHING
    amino_acid_df = amino_acid_df.sort_values(by=['codon_number']).rename(columns={'frequency': 'relative_frequency'})
    #FINDS GC PERCENTAGE OF EACH CODON
    amino_acid_df['gc_percentage'] = amino_acid_df['codon'].apply(GC)

    #RELATIVE FREQUENCY FINDS PROPORTION OF EACH CODON USE FOR EACH CODON POSITION
    total_counts = amino_acid_df.groupby('codon_number')['relative_frequency'].sum()
    amino_acid_df['relative_frequency'] = amino_acid_df['relative_frequency'] / amino_acid_df['codon_number'].map(total_counts)

    amino_acid_df.to_csv(output_filename, index=False)

fasta_file = '/Users/evehallett/Documents/dissertation/from_laurence/GCF_001043215.1_cds_from_genomic_refseq_filtered.fa'

# SERINE#
serine_dict = {"AGC", "AGT", 'TCT', 'TCC', 'TCA', 'TCG'}
process_amino_acid(fasta_file, serine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/serine.csv')

# LEUCINE#
leucine_dict = {'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'}
process_amino_acid(fasta_file, leucine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/leucine.csv')

#ARGININE#
arginine_dict = {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'}
process_amino_acid(fasta_file, arginine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/arginine.csv')

#ALANINE#
alanine_dict = {'GCT', 'GCC', 'GCA', 'GCG'}
process_amino_acid(fasta_file, alanine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/alanine.csv')

#GLYCINE#
glycine_dict = {'GGT', 'GGC', 'GGA', 'GGG'}
process_amino_acid(fasta_file, glycine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/glycine.csv')

#PROLINE#
proline_dict = {'CCT', 'CCC', 'CCA', 'CCG'}
process_amino_acid(fasta_file, proline_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/proline.csv')

#THREONINE#
threonine_dict = {'ACT', 'ACC', 'ACA', 'ACG'}
process_amino_acid(fasta_file, threonine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/threonine.csv')

#VALINE#
valine_dict = {'GTT', 'GTC', 'GTA', 'GTG'}
process_amino_acid(fasta_file, valine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/valine.csv')

#ISOLEUCINE#
isoleucine_dict = {'ATT', 'ATC', 'ATA'}
process_amino_acid(fasta_file, isoleucine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/isoleucine.csv')

#ASPARAGINE#
asparagine_dict = {'AAT', 'AAC'}
process_amino_acid(fasta_file, asparagine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/asparagine.csv')

#GLUTAMATE#
glutamate_dict = {'GAA', 'GAG'}
process_amino_acid(fasta_file, glutamate_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/glutamate.csv')

#ASPARTATE#
aspartate_dict = {'GAT', 'GAC'}
process_amino_acid(fasta_file, aspartate_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/aspartate.csv')

#CYSTINE#
cystine_dict = {'TGT', 'TGC'}
process_amino_acid(fasta_file, cystine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/cystine.csv')

#GLUTAMINE#
glutamine_dict = {'CAA', 'CAG'}
process_amino_acid(fasta_file, glutamine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/glutamine.csv')

#HISTIDINE#
histidine_dict = {'CAT', 'CAC'}
process_amino_acid(fasta_file, histidine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/histidine.csv')

#LYSINE#
lysine_dict = {'AAA', 'AAG'}
process_amino_acid(fasta_file, lysine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/lysine.csv')

#PHENYLALANINE#
phenylalanine_dict = {'TTT', 'TTC'}
process_amino_acid(fasta_file, phenylalanine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/phenylalanine.csv')

#TYROSINE#
tyrosine_dict = {'TAT', 'TAC'}
process_amino_acid(fasta_file, tyrosine_dict, '/Users/evehallett/Documents/dissertation/fig3/amino_acids_final/tyrosine.csv')

