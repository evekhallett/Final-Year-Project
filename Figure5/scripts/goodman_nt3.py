import pandas as pd
import os


#filter desired columns from Goodman (2013) data
data = pd.read_csv('/Users/evehallett/Documents/dissertation/to_submit/Figure5/goodman_data.csv', usecols=['Name', 'CDS.seq', 'Fltr.SetGood', 'Prot.FCC'])

#extract CDS
sequences = data['CDS.seq']
sequence_dfs_list = []

for seq in sequences:
    #codons 2 to 11
    codon_df_dict = {f'Codon_{i}': pd.DataFrame() for i in range(2, 12)}
    for base_position in range(3, 31, 3):
        if base_position + 2 < len(seq):
            current_row = data[data['CDS.seq'] == seq].copy()
            codon_position = (base_position + 3) // 3
            #finds the third nucleotide in each codon
            current_row[f'Position_{codon_position}'] = seq[base_position + 2].upper()
            codon_df_dict[f'Codon_{codon_position}'] = pd.concat([codon_df_dict[f'Codon_{codon_position}'],
                                                                current_row[['Name', 'CDS.seq', 'Fltr.SetGood', 'Prot.FCC', f'Position_{codon_position}']]], axis=1)


    sequence_df = pd.concat(codon_df_dict.values(), axis=1)
    sequence_df = sequence_df.loc[:, ~sequence_df.columns.duplicated()]
    sequence_dfs_list.append(sequence_df)

#adds filtered columns and the third site nucleotide for codons 2-11
result_df = pd.concat(sequence_dfs_list, ignore_index=True)
print(result_df.head())
os.chdir('/Users/evehallett/Documents/dissertation/to_submit/Figure5/files/')
result_df.to_csv('goodman_nt3.csv', index = False)
