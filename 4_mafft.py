import os
import subprocess

os.chdir("/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/aligned_sequences")
def align_sequences(input_file, output_file):
    # CHANGE MAFFT SETTING IF NEEDED
    subprocess.run(['/usr/local/bin/mafft', '--auto', input_file], stdout=open(output_file, 'w'))

def process_fasta_files(folder_path):
    # Iterate over each fasta file in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta"):
            file_path = os.path.join(folder_path, filename)

            # CHANGES FILE NAMES TO HAVE ALIGNED AT START
            output_file = os.path.join(folder_path, f"aligned_{filename}")
            align_sequences(file_path, output_file)

            print(f"Alignment completed for {filename}")

# Provide the path to the folder containing FASTA files
fasta_folder_path = "/Users/evehallett/Documents/dissertation/to_submit/Figure4/files/unaligned_sequences"

# Run the pipeline to align sequences in each FASTA file
process_fasta_files(fasta_folder_path)
