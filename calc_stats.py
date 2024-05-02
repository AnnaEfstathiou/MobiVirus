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
import os
from collections import Counter
from typing import Dict
from itertools import combinations
from Bio import SeqIO


"""
------------------
INTERNAL FUNCTIONS
------------------
""" 

def __csv_to_fasta(csv_file, output_fasta_file):

    """ Convert a CSV file to a FASTA file while excluding rows that are entirely NaN. """
    
    try:
        df = pd.read_csv(csv_file)  # Read the csv file
        with open(output_fasta_file, 'w') as fasta_file:
            for index, row in df.iterrows():
                if row.dropna().empty:
                    continue  # Skip rows that are all NaN
                # Convert each value to integer if it's a whole number, else keep original
                sequence = ''.join(str(int(val)) if isinstance(val, float) and val.is_integer() else str(val) for val in row)
                header = f'>{index + 1}'
                fasta_file.write(f'{header}\n{sequence}\n')
    except FileNotFoundError:
        print(f"Error: The file {csv_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

def __process_fasta(fasta_file, output_fasta):

    """ Processes a FASTA file to exclude rows that are entirely NaN. """

    try:
        # open input and output files
        with open(fasta_file, "r") as input_handle, open(output_fasta, "w") as output_handle:
            # parse the input FASTA file
            for line in SeqIO.parse(input_handle, "fasta"):
                # check if the sequence contains at least one position with 1
                if not 'nan' in line.seq:
                    # write the sequence to the output file
                    SeqIO.write(line, output_handle, "fasta")
    except FileNotFoundError:
        print(f"Error: The file {fasta_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

def __read_fasta(fasta_file):

    """ Function that creates a dictionary with headers as keys and sequences as values. """

    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                current_sequence = line.strip()[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()
    return sequences


def __check(sequences: Dict[str, str]) -> None:

    """ Check for valid sequence proportions. """
    """ Number of sequences > 1 """
    """ Sequences with the same length """

    if len(sequences) < 2:
        raise ValueError("At least 2 sequences required!")

    sequences_list = list(sequences.values()) # list of sequences
    seqs_len = len(sequences_list[0]) # sequence's length 
    if not all(len(s) == seqs_len for s in sequences_list):
        raise ValueError("Sequences are required to have the same length!")
    
def __polymorphic_sites(sequences: Dict[str, str]) -> int:
    
    """ Counts the number of positions which show differences (polymorphisms). """
    
    sequences_list = list(sequences.values())
    seg_sites = 0
    for i in range(len(sequences_list[0])):
        s = sequences_list[0][i]  # reference character
        if any(seq[i] != s for seq in sequences_list):
            seg_sites += 1
    return seg_sites

def __harmonic(n: int) -> float:

    """ Computes the n-1th harmonic number. """

    return sum(1 / i for i in range(1, n))

"""
----------
STATISTICS
----------
""" 

def tajimas_d(sequences: Dict[str, str]) -> float:

    """ Computes Tajima's D. """

    __check(sequences)

    seg_sites = __polymorphic_sites(sequences)
    if seg_sites == 0:
        return 0

    theta_pi = pi_estimator(sequences, safe=False)
    num_seq = len(sequences)
    harmonic = __harmonic(num_seq)
    a2 = sum(1 / (i**2) for i in range(1, num_seq))

    b1 = (num_seq + 1) / (3 * (num_seq - 1))
    b2 = (2 * (num_seq**2 + num_seq + 3)) / (9 * num_seq * (num_seq - 1))

    c1 = b1 - 1 / harmonic
    c2 = b2 - ((num_seq + 2) / (harmonic * num_seq)) + (a2 / (harmonic**2))

    if c1 == 0 or c2 == 0:
        return float('nan')  # informing that calculation is not reliable

    e1 = c1 / harmonic
    e2 = c2 / (harmonic**2 + a2)

    delta_Theta = theta_pi - (seg_sites / harmonic)
    tD = delta_Theta / ((e1 * seg_sites + e2 * seg_sites * (seg_sites - 1)) ** 0.5)
    return float(tD)

def pi_estimator(sequences: Dict[str, str], safe=True) -> float:

    """ Computes Pi estimatorn (Θπ). """

    if safe:
        __check(sequences)

    sequences_list = list(sequences.values())
    pairwise_combinations = combinations(sequences_list, 2)
    cs = [sum(not char1 == char2 for char1, char2 in zip(seq1, seq2)) for seq1, seq2 in pairwise_combinations]
    n = len(sequences_list)
    binomial = ((n - 1) * n) / 2  

    return sum(cs) / binomial


def watterson_estimator(sequences: Dict[str, str], safe=True) -> float:

    """ Computes Watterson estimator (Θw). """

    if safe:
        __check(sequences)

    seg_sites = __polymorphic_sites(sequences)
    harmonic = __harmonic(len(sequences))

    return seg_sites / harmonic


def count_haplotypes(sequences):

    """ Calculates the number of unique sequences. """

    return len(set(sequences.values()))


def calculate_haplotype_diversity(sequences):

    """ Calculates the diversity of haplotypes in a given set of sequences. """
    """ Formula: 1 - Σ(pi^2), pi is the frequency of the i-th haplotype divided by the total number of haplotypes. """

    total_haplotypes = len(sequences)
    haplotype_frequencies = Counter(sequences.values()) # count the frequencies of each unique haplotype 
    haplotype_diversity = 1 - sum((freq / total_haplotypes) ** 2 for freq in haplotype_frequencies.values()) # calculates the haplotype diversity based on the frequencies of each haplotype.
    return haplotype_diversity

"""
-------------
MAIN FUNCTION
-------------
""" 

def statistics(input_file):
    
    file_extension = os.path.splitext(input_file)[1] # get the file extension of the input file

    if file_extension.lower() == '.csv': # check if the input file is a CSV file

        new_fasta = input_file.rsplit('.', 1)[0] + '_processed.fa' # generate the processed FASTA file path by replacing the extension
        __csv_to_fasta(input_file, new_fasta) # convert CSV to FASTA (exclude nan sequences)

        sequences = __read_fasta(new_fasta) # create a dictionary with sequences 

        tajimas_d_score = tajimas_d(sequences)  # calculate Tajima's D
        pi_estimator_score = pi_estimator(sequences) # calculate Pi estimator
        watterson_estimator_score = watterson_estimator(sequences) # calculate Watterson estimator
        unique_count = count_haplotypes(sequences) # calculate the number of unique sequences
        haplotype_diversity = calculate_haplotype_diversity(sequences) # calculates the diversity of haplotypes
        
        os.remove(new_fasta) # delete the generated FASTA file

    elif file_extension.lower() == '.fasta' or file_extension.lower() == '.fa': # check if the input file is a FASTA file
   
        new_fasta = input_file.rsplit('.', 1)[0] + '_processed.fa' # generate the processed FASTA file path by replacing the extension
        __process_fasta(input_file, new_fasta) # exclude nan sequences from FASTA

        sequences = __read_fasta(new_fasta) # create a dictionary with sequences 

        tajimas_d_score = tajimas_d(sequences)  # calculate Tajima's D
        pi_estimator_score = pi_estimator(sequences) # calculate Pi estimator
        watterson_estimator_score = watterson_estimator(sequences) # calculate Watterson estimator
        unique_count = count_haplotypes(sequences) # calculate the number of unique sequences
        haplotype_diversity = calculate_haplotype_diversity(sequences) # calculates the diversity of haplotypes
        
        os.remove(new_fasta) # delete the generated FASTA file
    
    else:
        raise ValueError('The input file must be either a CSV or a FASTA file.')

    return tajimas_d_score, pi_estimator_score, watterson_estimator_score, unique_count, haplotype_diversity, len(sequences)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome CSV or FASTA file to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('input_file', type=str, help='The path to the input file (either CSV or FASTA).')
    args = parser.parse_args()
    input_file = args.input_file

    tajimas_d_score, pi_estimator_score, watterson_estimator_score, unique_count, haplotype_diversity, num_seqs  = statistics(input_file)
    
    print("Tajima's D score:", tajimas_d_score)
    print("Pi-Estimator score:", pi_estimator_score)
    print("Watterson-Estimator score:", watterson_estimator_score)
    print(f"Number of unique sequences: {unique_count}/{num_seqs} ({unique_count/num_seqs:.3f})")
    print("Haplotype diversity:", haplotype_diversity)


    