import re

def filter_fasta(input_file, output_file):
    #START CODONS = NTG
    #STOP CODONS = TAA, TGA, TAG, makes sure they are only at the end
    start_codon_pattern = re.compile(r'[A-Za-z][Tt][Gg]')
    stop_codon_pattern = re.compile(r'(TAA|TGA|TAG)[^ATCG]*$')
    #TO COUNT NO. OF SEQUENCES
    input_sequence_count = 0
    output_sequence_count = 0

    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        current_sequence = ''
        for line in in_file:
            if line.startswith('>'):  #SEQUENCE COUNT
                input_sequence_count += 1

                if (
                    len(current_sequence) % 3 == 0
                    and start_codon_pattern.search(current_sequence)
                    and stop_codon_pattern.search(current_sequence)
                    and all(base.upper() in 'ATCG' for base in current_sequence)
                ):
                    out_file.write(header + current_sequence + '\n')
                    output_sequence_count += 1

                header = line
                current_sequence = ''
            else:
                current_sequence += line.strip()
        #FOR LAST SEQUENCE
        if (
            len(current_sequence) % 3 == 0
            and start_codon_pattern.search(current_sequence)
            and stop_codon_pattern.search(current_sequence)
            and all(base.upper() in 'ATCG' for base in current_sequence)
        ):
            out_file.write(header + current_sequence + '\n')
            output_sequence_count += 1

    print(f'The input FASTA file has {input_sequence_count} sequences.')
    print(f'The output FASTA file has {output_sequence_count} sequences.')

#INPUT AND OUTPUT FASTAS
filter_fasta('/Users/evehallett/Documents/dissertation/from_laurence/GCF_001043215.1_cds_from_genomic_refseq.fa', '/Users/evehallett/Documents/dissertation/from_laurence/GCF_001043215.1_cds_from_genomic_refseq_filtered.fa')
