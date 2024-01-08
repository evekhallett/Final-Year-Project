import pandas as pd
import os

# READ IN CSV CONTAINING NUCLEOTIDE CONTENT AT THIRD SITE
os.chdir('/Users/evehallett/Documents/dissertation/to_submit/Figure7/files/')
data = pd.read_csv('goodman_nt3.csv')

#CODONS 2-11
base_columns = [f'Position_{i}' for i in range(2, 12)]
data[base_columns] = data[base_columns].applymap(str.upper)

# BINARY SCORING SYSTEM - GC = 1, AT = 0
binary_columns = data[base_columns].applymap(lambda base: 1 if base in ['G', 'C'] else 0)

# MAKE CSV
binary_data = pd.concat([data[['Name', 'Prot.FCC']], binary_columns], axis=1)
binary_data.to_csv('goodman_gc3.csv', index=False)
