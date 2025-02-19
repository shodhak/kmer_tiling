import pandas as pd
from Bio import SeqIO
import glob
import os
import time
import argparse

def read_fasta(fasta_file):
    # Read fasta file containing epitope sequences
    epitope_seq = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Save sequence names in a list
    seq_name = [record.id for record in epitope_seq]
    
    # Save sequences in a list
    sequences = [str(record.seq) for record in epitope_seq]
    
    # Save these in a data frame
    seq_df = pd.DataFrame({'annotation': seq_name, 'PROBE_SEQUENCE': sequences})
    
    # Add a column with epitope name
    seq_df['source'] = os.path.basename(fasta_file).split('.')[0]
    
    # Change all columns to string
    seq_df = seq_df.astype(str)
    
    # Return data frame
    return seq_df

def process_all_fastas(directory):
    start_time = time.time()
    
    # Find all fasta files in the specified directory
    fasta_files = glob.glob(os.path.join(directory, '*.fasta')) + \
                  glob.glob(os.path.join(directory, '*.faa')) + \
                  glob.glob(os.path.join(directory, '*.fa'))
    
    # Load each fasta file into a data frame and combine them
    combined_df = pd.concat([read_fasta(fasta_file) for fasta_file in fasta_files], ignore_index=True)
    
    end_time = time.time()
    print(f"Execution time for processing FASTA files: {end_time - start_time} seconds")
    
    return combined_df

def sliding_window(seq_string, probe_length, probe_offset):
    # Initialize empty list to store probes
    seq_probe = []
    for i in range(0, len(seq_string), probe_offset):
        if i + probe_length <= len(seq_string):
            tmp_probe = seq_string[i:i + probe_length]
            seq_probe.append(tmp_probe)
    return seq_probe

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA files and generate k-mers.")
    parser.add_argument("directory", type=str, help="Directory containing the FASTA files")
    parser.add_argument("probe_length", type=int, help="Length of the probes")
    parser.add_argument("probe_offset", type=int, help="Offset for the sliding window")
    args = parser.parse_args()

    # Directory containing the fasta files
    directory = args.directory

    # Process all fasta files and combine them into a single data frame
    combined_df = process_all_fastas(directory)

    # Apply sliding window function to each sequence in the combined data frame
    probe_length = args.probe_length
    probe_offset = args.probe_offset

    start_time = time.time()
    combined_df['PROBES'] = combined_df['PROBE_SEQUENCE'].apply(lambda x: sliding_window(x, probe_length, probe_offset))
    end_time = time.time()
    print(f"Execution time for generating probes: {end_time - start_time} seconds")

    # Save the sequences to a csv file
    start_time = time.time()
    combined_df[['annotation', 'PROBE_SEQUENCE', 'source']].to_csv(os.path.join(directory, 'combined_fasta_sequences.csv'), index=False)
    end_time = time.time()
    print(f"Execution time for saving sequences CSV: {end_time - start_time} seconds")

    # Explode the probes list into separate rows and save to a csv file
    start_time = time.time()
    probes_df = combined_df[['annotation', 'PROBES', 'source']].explode('PROBES')
    probes_df.to_csv(os.path.join(directory, 'combined_fasta_probes.csv'), index=False)
    end_time = time.time()
    print(f"Execution time for saving probes CSV: {end_time - start_time} seconds")