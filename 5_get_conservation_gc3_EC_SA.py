from Bio import SeqIO
import os

fasta_folder_path = "/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/aligned_sequences"
def calculate_percentage_difference(seq1, seq2):
    indels = sum(1 for base1, base2 in zip(seq1, seq2) if base1 == '-' or base2 == '-')
    if indels % 3 != 0:
        return 0
    differences = sum(1 for base1, base2 in zip(seq1, seq2) if (base1, base2) in [('G', 'A'), ('G', 'T'), ('C', 'A'), ('C', 'T'),('g', 'a'), ('g', 't'), ('c', 'a'), ('c', 't')])
    if differences >= 1:
        percentage_difference = differences * 100
    else:
        percentage_difference = 0
    return percentage_difference

def process_fasta_files(folder_path, codon_positions, min_seq_length = 120):
    result_dict = {pos: [] for pos in codon_positions}

    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta"):
            file_path = os.path.join(folder_path, filename)
            records = list(SeqIO.parse(file_path, "fasta"))

            if len(records) == 2:
                seq1 = records[0].seq #S. ENTERICA
                seq2 = records[1].seq #E. COLI

                for pos in codon_positions:
                    if len(seq1) % 3 == 0 and len(seq2) % 3 == 0 and len(seq1) >= min_seq_length and len(seq2) >= min_seq_length:
                        percentage_diff = calculate_percentage_difference(
                            seq1[(pos - 1) * 3 + 2],  # ACCESS 3RD NUCLEOTIDE SITE
                            seq2[(pos - 1) * 3 + 2]
                        )
                        result_dict[pos].append(percentage_diff)

    return result_dict

def write_output(result_dict, output_folder, output_filename):
    output_file_path = os.path.join(output_folder, output_filename)
    with open(output_file_path, 'w') as f:
        #COLUMS FOR THE CODON NUMBER AND PERECENTAGE CONSERVATION
        f.write("codon_number,percentage_conservation\n")
        for pos, values in result_dict.items():
            average_percentage_diff = sum(values) / len(values)
            #CALCULATES THE CONSERVATION %
            f.write(f"{pos},{100-average_percentage_diff}\n")

output_folder_path = "/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/"
output_filename = 'EC_salmonella_conservation_gc3.csv'

# CODON POSITION RANGE
codon_positions = list(range(2, 41))

# CALL FUNCTIONS
result_dict = process_fasta_files(fasta_folder_path, codon_positions)
write_output(result_dict, output_folder_path, output_filename)
