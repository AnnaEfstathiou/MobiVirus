import os
from typing import List, Dict
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse
from multiprocessing import Pool

from preprocessing import __csv_to_fasta, __read_fasta
from statistics import watterson_estimator, pi_estimator


def process_window_chunk(args):
    sequences, start, end, window_size = args
    theta_w_list = []
    pi_est_list = []

    for i in range(start, end):
        window_sequences = {key: seq[i:i+window_size] for key, seq in sequences.items()}
        theta_w = watterson_estimator(window_sequences)  # calculate Watterson estimator
        theta_w_list.append(theta_w)                     # append the score to a list
        pi_est = pi_estimator(window_sequences)          # calculate Pi estimator
        pi_est_list.append(pi_est)                       # append the score to a list

    return theta_w_list, pi_est_list


def sliding_window_watterson(sequences: Dict[str, str], window_size: int, num_processes: int = 8) -> List[float]:
    if not sequences:
        raise ValueError("No sequences found. Please check the input FASTA file.")
    
    sequence_length = len(next(iter(sequences.values())))
    num_windows = sequence_length - window_size + 1
    chunk_size = num_windows // num_processes

    args = [(sequences, i, min(i + chunk_size, num_windows), window_size) for i in range(0, num_windows, chunk_size)]

    with Pool(num_processes) as pool:
        results = pool.map(process_window_chunk, args)

    theta_w_list = [item for sublist in [result[0] for result in results] for item in sublist]
    pi_est_list = [item for sublist in [result[1] for result in results] for item in sublist]

    return theta_w_list, pi_est_list, sequence_length


def plot_selective_sweep(theta_w: List[float], pi_est: List[float], window_size: int, num_sequences, sequence_length):
    plt.figure(figsize=(10, 6))
    plt.plot(theta_w, marker='.', linestyle='solid', color='navy', label='Watterson Estimator')
    plt.plot(pi_est, marker='.', linestyle='solid', color='firebrick', label='Pi Estimator')
    plt.axvline(x=4999, color='green', linestyle='dashed', label='SNP position')  # add a green line at SNP position
    plt.plot([], [], ' ', label=f'Number of Sequences: {num_sequences}')
    plt.plot([], [], ' ', label=f'Sequence Length: {sequence_length}')
    plt.xlabel('Window Position')
    plt.ylabel('Estimator Value')
    plt.title(f'Selective Sweep with Window Size {window_size}')
    plt.legend()
    if args.save_png:
        plt.savefig(args.save_png, format="png")
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
            
        theta_w, pi_est, seq_len = sliding_window_watterson(sequences, window_size) # Create 2 list with the corresponding scores
        plot_selective_sweep(theta_w, pi_est, window_size, len(sequences), seq_len) # Plot scores
        
        os.remove(output_fasta_file) # Delete the temporary created fasta file

    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome CSV or FASTA file to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('genome_file', type=str, help='The path to the input CSV file.')
    parser.add_argument('-w','--window_size', type=int, help='Length of sliding window.', default=None)
    parser.add_argument('-s','--save_png', type=str, help='Flag to save the plot as an PNG file.')
    args = parser.parse_args()

    main(args.genome_file, args.window_size)
