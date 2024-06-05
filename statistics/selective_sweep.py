""" 
INPUT: CSV file with binary sequences
Calculating Pi-Estimator score & Watterson-Estimator score for sliding windows.
Plot the scores to depict the selective sweep.
"""

import os
from typing import List, Dict
import matplotlib.pyplot as plt
import argparse

from preprocessing import __csv_to_fasta, __read_fasta
from statistics import watterson_estimator, pi_estimator


def sliding_window_watterson(sequences: Dict[str, str], window_size: int) -> List[float]:

    """ Calculate the Watterson estimator and Pi estimator for each sliding window. """

    if not sequences:
        raise ValueError("No sequences found. Please check the input FASTA file.")
    
    sequence_length = len(next(iter(sequences.values())))
    theta_w_list = []
    pi_est_list = []

    for i in range(sequence_length - window_size + 1):
        window_sequences = {key: seq[i:i+window_size] for key, seq in sequences.items()}

        theta_w = watterson_estimator(window_sequences) # calculate pi estimator
        theta_w_list.append(theta_w)                    # append the score to a list

        pi_est = pi_estimator(window_sequences)         # calculate watterson estimator
        pi_est_list.append(pi_est)                      # append the score to a list

    return theta_w_list, pi_est_list

def plot_watterson(theta_w: List[float], pi_est: List[float], window_size: int):

    """ Plot the Watterson estimator and Pi estimator values. """

    plt.figure(figsize=(10, 6))
    plt.plot(theta_w, marker='o', linestyle='--', color='b', label='Watterson Estimator')
    plt.plot(pi_est, marker='o', linestyle='--', color='r', label='Pi Estimator')
    plt.title(f'Selective Sweep (Window Size: {window_size})')
    plt.xlabel('Window Position')
    plt.ylabel('Scores')
    plt.grid(True)
    plt.legend()
    if args.save_svg:
        plt.savefig("selective_sweep.svg", format="svg")
    else:
        plt.show()
    plt.close()

def main(csv_file, window_size):

    output_fasta_file = 'fasta_file.fa'
    __csv_to_fasta(csv_file, output_fasta_file) 

    try:
        sequences = __read_fasta(output_fasta_file) 
        if not sequences:
            print("No sequences loaded. Please check the FASTA file.")
        else:
            print(f"Loaded {len(sequences)} sequences.")
            
        theta_w, pi_est = sliding_window_watterson(sequences, window_size) # Create 2 list with the corresponding scores
        plot_watterson(theta_w, pi_est, window_size) # Plot scores
        
        os.remove(output_fasta_file) # Delete the temporary created fasta file

    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome CSV or FASTA file to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('genome_file', type=str, help='The path to the input CSV file.')
    parser.add_argument('-w','--window_size', type=int, help='Length of sliding window.', default=None)
    parser.add_argument('-s','--save_svg', action='store_true', help='Flag to save the plot as an SVG file.')
    args = parser.parse_args()

    main(args.genome_file, args.window_size)