""" 
INPUT: FASTA or CSV file with binary sequences
Calculating: 
- Tajima's D 
- Pi-Estimator score
- Watterson-Estimator score
- number of unique sequences
- haplotype diversity 
"""

import argparse
import pandas as pd
from collections import Counter
import subprocess
from Bio import SeqIO
import os
from _tajimas_d import *


# def csv_to_fasta_and_preprocess(csv_file, output_fasta_file):

#     """ Convert a CSV file directly to a preprocessed FASTA file, excluding all-zero sequences. """
    
#     df = pd.read_csv(csv_file) # read the csv file as 
#     with open(output_fasta_file, 'w') as fasta_file:
#         for index, row in df.iterrows():
#             sequence = ''.join([str(int(val)) if val.is_integer() else str(val) for val in row.values])
#             if '1' in sequence:  # check if the sequence contains at least one position with 1
#                 header = f'>{index + 1}'
#                 fasta_file.write(f'{header}\n{sequence}\n')

def csv_to_fasta_and_preprocess(csv_file, output_fasta_file):

    """ Convert a CSV file directly to a preprocessed FASTA file, excluding all-zero sequences. """
    
    df = pd.read_csv(csv_file) # read the csv file as 
    with open(output_fasta_file, 'w') as fasta_file:
        for index, row in df.iterrows():
            sequence = ''.join([str(int(val)) if val.is_integer() else str(val) for val in row.values])
            # if '1' in sequence:  # check if the sequence contains at least one position with 1
            #     header = f'>{index + 1}'
            #     fasta_file.write(f'{header}\n{sequence}\n')
            header = f'>{index + 1}'
            fasta_file.write(f'{header}\n{sequence}\n')


def preprocess_fasta(input_fasta, output_fasta):

    """ Preprocesses the FASTA file to remove lines where all positions are 0. """

    # open input and output files
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        # parse the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # check if the sequence contains at least one position with 1
            if '1' in record.seq:
                # write the sequence to the output file
                SeqIO.write(record, output_handle, "fasta")


def read_fasta(file_path):

    """ Function that creates a dictionary with headers as keys and sequences as values"""

    sequences = {}
    current_sequence = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                current_sequence = line.strip()[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()
    return sequences


def run_tajimas_d(fasta_file_path):

    """ Function to run Tajima's D on a FASTA file """

    safe_fasta_file_path = f"'{fasta_file_path}'" # ensure the file path is safely quoted to prevent shell interpretation issues
    tajimas_D_command = f'tajimas_d -f {safe_fasta_file_path} -p -t -w' # construct the command to run Tajima's D
    result = subprocess.run(tajimas_D_command, shell=True, capture_output=True, text=True) # execute the command in the shell and capture output
    # parsing the results to create a dict
    parsed_result = {}
    lines = result.stdout.split('\n')
    for line in lines:
        if line.strip(): # ignoring empty lines
            key, value = line.split(':')
            parsed_result[key.strip()] = value.strip()
    return parsed_result


def count_haplotypes(sequences):

    """ Function that calculates the number of unique sequences """

    return len(set(sequences.values()))


def calculate_haplotype_diversity(sequences):

    """ Function that calculates the diversity of haplotypes in a given set of sequences """
    """ Formula: 1 - Σ(pi^2), pi is the frequency of the i-th haplotype divided by the total number of haplotypes """

    total_haplotypes = len(sequences)
    haplotype_frequencies = Counter(sequences.values()) # count the frequencies of each unique haplotype 
    haplotype_diversity = 1 - sum((freq / total_haplotypes) ** 2 for freq in haplotype_frequencies.values()) # calculates the haplotype diversity based on the frequencies of each haplotype.
    return haplotype_diversity


def run_statistics_for_file(file_path):

    file_extension = os.path.splitext(file_path)[1] # get the file extension of the input file

    if file_extension.lower() == '.csv': # check if the input file is a CSV file

        preprocessed_file_path = file_path.rsplit('.', 1)[0] + '_processed.fa' # generate the processed FASTA file path by replacing the extension
        csv_to_fasta_and_preprocess(file_path, preprocessed_file_path) # remove lines where all positions are 0

        statistical_results = run_tajimas_d(preprocessed_file_path) # Tajima's D, Pi-Estimator score, Watterson-Estimator score
        sequences = read_fasta(preprocessed_file_path)
        unique_count = count_haplotypes(sequences) # calculate the number of unique sequences
        haplotype_diversity = calculate_haplotype_diversity(sequences) # calculates the diversity of haplotypes
        
        os.remove(preprocessed_file_path) # delete the generated FASTA file

    elif file_extension.lower() == '.fasta' or file_extension.lower() == '.fa': # check if the input file is a FASTA file
   
        preprocessed_file_path = file_path.rsplit('.', 1)[0] + '_processed.fa' # generate the processed FASTA file path by replacing the extension
        preprocess_fasta(file_path, preprocessed_file_path) # remove lines where all positions are 0

        statistical_results = run_tajimas_d(preprocessed_file_path) # Tajima's D, Pi-Estimator score, Watterson-Estimator score
        sequences = read_fasta(preprocessed_file_path)
        unique_count = count_haplotypes(sequences) # calculate the number of unique sequences
        haplotype_diversity = calculate_haplotype_diversity(sequences) # calculates the diversity of haplotypes

        os.remove(preprocessed_file_path) # delete the generated FASTA file
    
    else:
        raise ValueError('The input file must be either a CSV or a FASTA file.')

    return statistical_results, unique_count, haplotype_diversity, len(sequences)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome CSV or FASTA file to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('input_file', type=str, help='The path to the input file (either CSV or FASTA).')
    args = parser.parse_args()
    input_file = args.input_file

    statistical_results, unique_count, haplotype_diversity, num_seqs = run_statistics_for_file(input_file)  
    print(f"Number of unique sequences: {unique_count}/{num_seqs} ({unique_count/num_seqs:.4f})")
    print("Haplotype diversity:", haplotype_diversity)
    print("Pi-Estimator score:", statistical_results["Pi-Estimator score"])
    print("Watterson-Estimator score:", statistical_results["Watterson-Estimator score"])
    print("Tajima's D score:", statistical_results["Tajima's D score"])

