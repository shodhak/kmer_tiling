import pandas as pd
from Bio import SeqIO
import glob
import os
import time
import argparse
import bodo
import numpy as np


df_sample = pd.DataFrame(
    {"annotation": ["a"], "SEQUENCE": ["p"], "source": ["s"]}
)
df_type = bodo.typeof(df_sample)

@bodo.wrap_python(df_type)
def read_fasta(fasta_file):
    # Save sequence names and sequences in a list
    seq_names = []
    sequences = []

    # Read fasta file containing epitope sequences
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_names.append(record.id)
        sequences.append(str(record.seq))

    # Save these in a data frame
    seq_df = pd.DataFrame({"annotation": seq_names, "SEQUENCE": sequences})

    # Add a column with epitope name
    seq_df["source"] = os.path.basename(fasta_file).split(".")[0]

    # Change all columns to string
    seq_df = seq_df.astype(str)
    return seq_df

#@bodo.jit(cache=True)
def sliding_window(seq_string, probe_length, probe_offset):
    """Parallel sliding window function"""
    probes = []
    for i in range(0, len(seq_string), probe_offset):
        if i + probe_length <= len(seq_string):
            tmp_probe = seq_string[i:i + probe_length]
            probes.append(tmp_probe)
    return probes

#@bodo.jit(cache=True)
def process_sequences_parallel(sequences, probe_length, probe_offset):
    """Parallel processing of sequences to generate probes"""
    return [sliding_window(seq, probe_length, probe_offset) for seq in sequences]

#@bodo.jit(cache=True)
def process_all_fastas(directory):
    start_time = time.time()
    
    # Find all fasta files in the specified directory
    fasta_files = glob.glob(os.path.join(directory, '*.fasta')) + \
                  glob.glob(os.path.join(directory, '*.faa')) + \
                  glob.glob(os.path.join(directory, '*.fa'))
    
    # First, read all files sequentially (BioPython operations)
    dfs = []
    for fasta_file in fasta_files:
        df = read_fasta(fasta_file)
        dfs.append(df)
    
    combined_df = pd.concat(dfs, ignore_index=True)
    
    end_time = time.time()
    print(f"Execution time for processing FASTA files: {end_time - start_time} seconds")
    
    return combined_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA files and generate k-mers.")
    parser.add_argument("directory", type=str, help="Directory containing the FASTA files")
    parser.add_argument("probe_length", type=int, help="Length of the probes")
    parser.add_argument("probe_offset", type=int, help="Offset for the sliding window")
    parser.add_argument("--threads", type=int, default=None, 
                       help="Number of threads to use (default: use all available)")
    args = parser.parse_args()

    if args.threads:
        os.environ["BODO_NUM_THREADS"] = str(args.threads)

    # Process all fasta files and combine them into a single data frame
    combined_df = process_all_fastas(args.directory)

    # Generate probes in parallel
    start_time = time.time()
    sequences = combined_df['SEQUENCE'].tolist()
    probes_lists = process_sequences_parallel(sequences, args.probe_length, args.probe_offset)
    combined_df['PROBES'] = probes_lists
    end_time = time.time()
    print(f"Execution time for generating probes: {end_time - start_time} seconds")

    # Save the sequences to a csv file
    start_time = time.time()
    combined_df[['annotation', 'SEQUENCE', 'source']].to_csv(
        os.path.join(args.directory, 'combined_fasta_sequences.csv'), 
        index=False
    )
    end_time = time.time()
    print(f"Execution time for saving sequences CSV: {end_time - start_time} seconds")

    # Explode the probes list into separate rows and save to a csv file
    start_time = time.time()
    probes_df = combined_df[['annotation', 'PROBES', 'source']].explode('PROBES')
    probes_df.to_csv(os.path.join(args.directory, 'combined_fasta_probes.csv'), index=False)
    end_time = time.time()
    print(f"Execution time for saving probes CSV: {end_time - start_time} seconds")