from Bio import SeqIO
import pandas as pd
import numpy as np

def calculate_gc3_content(sequence):
    gc3_content = []
    #40 CODONS
    if len(sequence) >= 120:
        for i in range(5, 121, 3):
            if i + 2 < len(sequence):
                #THIRD POSITION IN CODON
                gc3 = sequence[i+2]
                if gc3.upper() in {'G', 'C'}:
                    gc3_value = 1
                else:
                    gc3_value = 0
                gc3_content.append(gc3_value)
            else:
                break
        return gc3_content
    else:
        return None

def main():
    file_path = '/Users/evehallett/Documents/dissertation/from_laurence/GCF_001043215.1_cds_from_genomic_refseq_filtered.fa'
    data = []
    #READ IN FASTA FILE
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence = str(record.seq)
            overall_gc3 = calculate_gc3_content(sequence)
            if overall_gc3 is not None:
                data.append(overall_gc3)
    #MAKING COLUMN NAMES CORRESPOND TO CODON POSITION
    columns = [f'codon_{i}' for i in range(2, 41)]
    df = pd.DataFrame(data, columns=columns)
    #DEALING WITH MISSING VALUES
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna()
    df = df.astype(int)
    #MAKE CSV
    df.to_csv('/Users/evehallett/Documents/dissertation/to_submit/Figure2/ecoli_gc3_40.csv', index=False)
    print(df)

if __name__ == "__main__":
    main()
